# NumDiff Code Review

**Scope:** all of `https://github.com/jacobwilliams/NumDiff/blob/master/src/` (≈6,100 lines across 6 modules), `tests/`, `fpm.toml`, and CI, at commit `954a785`.
**Method:** I read every source file. I then compiled the library with `gfortran 15.3 -g -fcheck=all` (pixi environment) and wrote a small reproducer for each high-severity finding. Every finding marked **✅ reproduced** comes with output from an actual run.

---

## Summary

The library is well documented, uses the FORD conventions consistently, and the ported algorithms (Oliver/Kahaner `diff`, Coleman–Garbow–Moré `dsm`) look faithful to the originals. The problems are in the newer glue code in `numerical_differentiation_module`:

| # | Severity | Issue | Status |
|---|----------|-------|--------|
| 1 | 🔴 High | Partitioned Jacobian is **silently wrong** when a linear sparsity pattern is used | ✅ reproduced |
| 2 | 🔴 High | Function cache returns **wrong Jacobian values** on partial cache hits | ✅ reproduced |
| 3 | 🔴 High | `class=11/13/15/17` crashes because the methods with IDs 500–800 are never found | ✅ reproduced |
| 4 | 🔴 High | `sparsity_mode=3` without calling `set_sparsity_pattern` causes a segfault | ✅ reproduced |
| 5 | 🟠 Medium | `num_sparsity_points` ≳ 164 reads out of bounds | ✅ reproduced |
| 6 | 🟠 Medium | `terminate()` does not stop `diff` mode: it does about 5× more work instead | ✅ reproduced |
| 7 | 🟠 Medium | Re-initializing an object does not clear exceptions or the old sparsity pattern | ✅ reproduced |
| 8 | 🟡 Low | `function_cache%print` has an off-by-one error (0-based table) | ✅ reproduced |
| — | 🟡 Low | Assorted input validation, numerical-robustness, and performance items | see below |

The tests contain **no assertions**, so none of the issues above would be caught by CI (see [Testing & CI](#testing--ci)).

---

## 🔴 High severity

### 1. Partitioned Jacobian is wrong when linear elements exist

**Where:** [numerical_differentiation_module.f90:1734-1753](https://github.com/jacobwilliams/NumDiff/blob/master/src/numerical_differentiation_module.f90#L1734-L1753) (and [1978](https://github.com/jacobwilliams/NumDiff/blob/master/src/numerical_differentiation_module.f90#L1978), [2172](https://github.com/jacobwilliams/NumDiff/blob/master/src/numerical_differentiation_module.f90#L2172)); computation in `compute_jacobian_partitioned` ([2753](https://github.com/jacobwilliams/NumDiff/blob/master/src/numerical_differentiation_module.f90#L2753)).

The DSM partition is computed from the **nonlinear** pattern (`irow`/`icol`) only. Linear (constant) elements are left out of that pattern, but the function still depends on those variables. So when two columns share a row through a linear element, DSM can put them in the same group. Perturbing both columns together then adds the linear term into the nonlinear element.

**Repro:** `f1 = x1² + 3·x2`, `f2 = x2·x3`. Nonlinear pattern `(1,1),(2,2),(2,3)`, linear `(1,2)=3`, forward difference, `partition_sparsity_pattern=.true.`:

```
row1:    5.000010    3.000000    0.000000     <-- J(1,1) should be 2
row2:    0.000000    3.000000    2.000000
```

Columns 1 and 2 share no row in the nonlinear pattern, so they end up in one group, and `∂f1/∂x1` picks up the `3·dx2/dx1` term.

**Fix:** Partition on the union of the nonlinear and linear patterns (pass the concatenated `irow`/`icol` to `dsm`), and keep using only the nonlinear indices when extracting `jac`.

> Note: `compute_linear_sparsity_pattern` is private, defaults to `.false.`, and **nothing ever sets it** (grep confirms). So the automatic linear detection in `compute_sparsity_random`/`_random_2` is currently dead code. Only `set_sparsity_pattern(..., linear_*)` reaches this bug today. If the flag is meant to be user-facing, add it to `initialize`; otherwise remove the dead branches.

### 2. Cache partial hits corrupt the returned `f`

**Where:** `compute_function_with_cache`, [numerical_differentiation_module.f90:389-397](https://github.com/jacobwilliams/NumDiff/blob/master/src/numerical_differentiation_module.f90#L389-L397)

On a partial hit, `cache%get` puts the cached elements into `f`. Then `problem_func(x, f, missing)` is called with that same `f`. The dummy argument is `intent(out)`, so the cached values become **undefined**. In practice they are lost as soon as the user's function initializes `f` (for example `f = 0`, which is a very common pattern).

**Repro:** the same function as #1, with `sparsity_mode=2`, `cache_size=100`, and a user function that starts with `f = 0`:

```
with cache row1:    2.000010************    0.000000     <-- J(1,2) ≈ 3e5, should be 3
```

Column 1 caches `f1(x)`. Column 2 asks for `[f1,f2](x)`, gets a partial hit, and `f1(x)` is wiped.

**Fix:** evaluate into a temporary and merge:

```fortran
real(wp),dimension(me%m) :: ftmp
...
call me%problem_func(x, ftmp, missing)
f(missing) = ftmp(missing)
call me%cache%put(i, x, ftmp, missing)
```

### 3. Classes 11, 13, 15, 17 are unreachable and crash

**Where:** `get_all_methods_in_class`, [numerical_differentiation_module.f90:809-835](https://github.com/jacobwilliams/NumDiff/blob/master/src/numerical_differentiation_module.f90#L809-L835); IDs defined at [753-777](https://github.com/jacobwilliams/NumDiff/blob/master/src/numerical_differentiation_module.f90#L753-L777).

The loop increments `id` from 1 and exits at the first ID that is not found (45). The 11-, 13-, 15- and 17-point methods have IDs 500, 600, 700 and 800, so they are never returned. `initialize(class=11)` succeeds, and `compute_jacobian` then fails:

```
init failed? F
At line 876 ...: Fortran runtime error: Allocatable argument 'list_of_methods' is not allocated
```

Any invalid class (for example 1 or 10) crashes the same way.

**Fix:**
- Iterate over an explicit list of valid IDs (for example a module-level `parameter` array), or make the IDs contiguous.
- In `initialize_numdiff`, raise an exception when `get_all_methods_in_class` returns an unallocated list.
- Guard `select_finite_diff_method*` against an unallocated `meth`.
- Update the doc comment at [512](https://github.com/jacobwilliams/NumDiff/blob/master/src/numerical_differentiation_module.f90#L512) ("id codes … are sequential"), which is no longer true.

### 4. `sparsity_mode=3` without a user pattern segfaults

**Where:** `compute_jacobian`, [numerical_differentiation_module.f90:2430](https://github.com/jacobwilliams/NumDiff/blob/master/src/numerical_differentiation_module.f90#L2430)

Mode 3 sets `compute_sparsity => null()`. If the user calls `compute_jacobian` before `set_sparsity_pattern`, the null procedure pointer is called and the program gets a SIGSEGV (reproduced).

**Fix:**
```fortran
if (.not. me%sparsity%sparsity_computed) then
    if (.not. associated(me%compute_sparsity)) then
        call me%raise_exception(30,'compute_jacobian','sparsity pattern has not been set.')
        return
    end if
    call me%compute_sparsity(x)
    if (me%exception_raised) return
end if
```
(This also adds the missing exception check after `compute_sparsity`.)

---

## 🟠 Medium severity

### 5. `divide_interval` returns fewer points than requested

**Where:** [utilities_module.f90:469](https://github.com/jacobwilliams/NumDiff/blob/master/src/utilities_module.f90#L469), used at [numerical_differentiation_module.f90:2088-2114](https://github.com/jacobwilliams/NumDiff/blob/master/src/numerical_differentiation_module.f90#L2088-L2114)

`delta*i*noise` with `noise ≈ 1.0123` exceeds 1 for `i ≳ 0.988·(n+1)`:
- For `n ≥ ~81`, the last point is clamped to the upper bound.
- For `n ≳ 164`, several points clamp to the same value, and `unique()` removes the duplicates.

`compute_sparsity_random_2` then indexes `coeffs(1:num_sparsity_points)`:

```
Fortran runtime error: Index '200' of dimension 1 of array 'coeffs' above upper bound of 199
```

**Fix:** scale so that the largest point stays interior, for example `tmp(i) = delta*i*(1 + (noise-1)*(1-2*delta))`, or apply the noise as `delta*(i + small_offset)`. Then assert `size(points)==num_points`. Also validate `num_sparsity_points >= 1` in `initialize`.

### 6. `terminate()` inside `diff` mode does not terminate

**Where:** `dfunc` in `compute_jacobian_with_diff`, [numerical_differentiation_module.f90:2735-2737](https://github.com/jacobwilliams/NumDiff/blob/master/src/numerical_differentiation_module.f90#L2735-L2737)

When an exception is raised, `dfunc` calls `me%terminate()` (the `numdiff_type` version, which is a no-op because an exception is already set) instead of `this%terminate()` (the `diff_func` version, which sets `ifail=-1`). `diff` keeps running on `fx = 0` values, and the outer loop moves on to every remaining element.

| scenario | user function calls |
|---|---|
| normal run | 804 |
| user calls `terminate()` at call 5 | **3850** |

**Fix:** `call this%terminate()`. The `ifail == -1` branch at [2699](https://github.com/jacobwilliams/NumDiff/blob/master/src/numerical_differentiation_module.f90#L2699) then works as its comment describes.

### 7. Re-initialization leaves stale state

**Where:** `initialize_numdiff` / `initialize_numdiff_for_diff`; `clear_exceptions` ([3160](https://github.com/jacobwilliams/NumDiff/blob/master/src/numerical_differentiation_module.f90#L3160)) is private and **never called**.

- **Exceptions are sticky.** After any exception (including a user `terminate()`), calling `initialize` again does nothing useful: every setter returns early on `exception_raised`. Reproduced: `failed after re-initialize: T`, and `jac` is unallocated. The only way out is `destroy()`, which the docs don't mention.
- **The sparsity pattern is not reset.** `sparsity%sparsity_computed` stays `.true.` across re-initialization, so a new `n`/`m`/`sparsity_mode` silently reuses the old pattern. `tests/test1.f90:235` already hits this: `diff_initialize(..., sparsity_mode=1)` never computes the dense pattern and reuses the partitioned pattern from the previous class-9 run. It only works because the function is the same.
- **Other fields persist too.** `partition_sparsity_pattern`, `print_messages`, `eps`/`acc`, `info_function`, `meth`/`class`, etc.

**Fix:** start both initializers with `call me%destroy()` (or at least `clear_exceptions` + `destroy_sparsity_pattern`). Take care: `destroy` is `intent(out)` on `me`, so the `problem_func` argument must not alias anything inside `me`. Also consider making `clear_exceptions` public.

Related: in the `classes` branch, the exception at [1450-1455](https://github.com/jacobwilliams/NumDiff/blob/master/src/numerical_differentiation_module.f90#L1450-L1455) is missing the `return` that the three sibling branches have.

---

## 🟡 Low severity / correctness hygiene

### Cache module (`cache_module.f90`)
- **Off-by-one in `print_cache`** ([98](https://github.com/jacobwilliams/NumDiff/blob/master/src/cache_module.f90#L98)): the table is allocated `0:isize-1`, but the loop runs `1..size(me%c)`. It skips slot 0 and reads past the end (reproduced: `Index '4' … above upper bound of 3`). Use `lbound/ubound`.
- **Wrong bounds check in `put_in_cache`** ([188](https://github.com/jacobwilliams/NumDiff/blob/master/src/cache_module.f90#L188)): `i<=size(me%c)` should be `i>=0 .and. i<size(me%c)`.
- **Hash is not portable** ([266](https://github.com/jacobwilliams/NumDiff/blob/master/src/cache_module.f90#L266)):
  - `ishft(hash,5)+hash+…` relies on signed-integer overflow, which the standard does not define. It traps under `-ftrapv`/sanitizers.
  - `abs(-huge-1)` also overflows.
  - For `REAL32`, `transfer(r(i), 1_ip)` produces an int64 whose upper 4 bytes are processor-dependent, so the hash can be non-deterministic and the cache misses.
  - For `REAL128`, only half the bits are hashed.
  - The comment at [8](https://github.com/jacobwilliams/NumDiff/blob/master/src/cache_module.f90#L8) ("same number of bits as real(wp)") is true only for `REAL64`.
  - Suggested fix: hash the bytes via `transfer(r, [0_int8])`, use FNV-1a with `iand`/masking, or compute modulo a prime at each step.
- `+0.0` and `-0.0` compare equal but hash differently. This is harmless (just a cache miss) but worth a comment.

### Utilities (`utilities_module.f90`)
- `unique_*` ([159](https://github.com/jacobwilliams/NumDiff/blob/master/src/utilities_module.f90#L159), [196](https://github.com/jacobwilliams/NumDiff/blob/master/src/utilities_module.f90#L196)) and `equal_within_tol` ([429](https://github.com/jacobwilliams/NumDiff/blob/master/src/utilities_module.f90#L429)) index element 1 without checking for an empty input.
- `expand_vector(..., finished=.true.)` with `vec` unallocated does `tmp = vec(1:n)` on an unallocated array ([79](https://github.com/jacobwilliams/NumDiff/blob/master/src/utilities_module.f90#L79)). Callers currently guard this, but the routine should too.
- `chunk_size=0` is accepted (`abs(chunk_size)` at [numerical_differentiation_module.f90:1091](https://github.com/jacobwilliams/NumDiff/blob/master/src/numerical_differentiation_module.f90#L1091), [1484](https://github.com/jacobwilliams/NumDiff/blob/master/src/numerical_differentiation_module.f90#L1484)). It makes `expand_vector` allocate size 0 and then write `vec(1)`. Use `max(1, abs(chunk_size))`.
- The module has no module-level `implicit none` (each procedure has its own). This is harmless, but inconsistent with the other modules.
- The `swap_real` doc comment says "Swap two integer values".

### Input validation (`numerical_differentiation_module.f90`)
- `set_sparsity_pattern` ([1722](https://github.com/jacobwilliams/NumDiff/blob/master/src/numerical_differentiation_module.f90#L1722)) checks only upper bounds. Zero or negative `irow`/`icol` pass and later index out of bounds. The same applies to the linear arrays. Duplicate pairs are not rejected either.
- A user-supplied `ngrp`/`maxgrp` is range-checked but not checked for **consistency** (two columns in one group sharing a row). A cheap O(nnz) check would prevent silently wrong Jacobians.
- `dpert` is never checked for `> 0`. With `dpert=0`, `perturb_mode=1` gives `dx=0` and a division by zero in `df/(den*dx)`. The fallback `where (dx<eps) dx = dpert` at [2959](https://github.com/jacobwilliams/NumDiff/blob/master/src/numerical_differentiation_module.f90#L2959) doesn't help when `dpert` itself is tiny or zero.
- `n`, `m`, `num_sparsity_points`, and `cache_size` are not validated for positivity.

### `diff_module.f90`
- `faccur` is called with `f1` still undefined ([281](https://github.com/jacobwilliams/NumDiff/blob/master/src/diff_module.f90#L281)). This is benign because `h0=0` takes the branch that ignores `f1`. Initialize `f1` anyway to keep `-Wmaybe-uninitialized`/`-finit-real=snan` runs clean.
- `deriv`/`error` are undefined on `ifail = -1, 2, 3`. Consider setting them (for example to `0` and `huge`) so callers can't pick up garbage.

### `dsm_module.f90`
- `iwa(max(m,6*n))` ([97](https://github.com/jacobwilliams/NumDiff/blob/master/src/dsm_module.f90#L97)) is an automatic (stack) array. For large sparse problems (for example `n` ≈ 10⁶ gives 24 MB) this will overflow the default stack, the same class of problem as PR #12. Make it `allocatable`. **FIXED**
- In `ido`, `maxlst = maxlst/n` ([398](https://github.com/jacobwilliams/NumDiff/blob/master/src/dsm_module.f90#L398)) can be 0 for very sparse matrices (Σ row_nnz² < n). The selection loop then runs zero times and `jcol` is stale. This is inherited from MINPACK and I couldn't construct a case that reaches `ido` with that pattern, but `max(1, …)` is a free safety net.
- Many helper arguments (`seq`, `slo`, `numsrt`, `srtdat`, `fdjs`) lack `intent`. Adding intents would let the compiler catch misuse.

---

## Numerical / design observations

1. **Absolute tolerances for sparsity detection.** `function_precision_tol` and `linear_sparsity_tol` both default to `epsilon(1.0_wp)` and are compared **absolutely** ([1943](https://github.com/jacobwilliams/NumDiff/blob/master/src/numerical_differentiation_module.f90#L1943), [2142](https://github.com/jacobwilliams/NumDiff/blob/master/src/numerical_differentiation_module.f90#L2142)).
   - For functions of magnitude around 1e6, round-off noise alone marks elements as nonzero.
   - For functions of magnitude around 1e-20, real dependencies are dropped.

   A relative test (`abs(a-b) <= tol*max(1,abs(a),abs(b))`) is more robust. At minimum, document that users with badly scaled functions must set these tolerances.
2. **Step representability.** `xp(column) = x + dx_factor*dx` is used, but the divisor is `dx`, not the step that was actually representable, `(x+h) - x`. For large `|x|` this loses accuracy, and the standard trick is cheap.
3. **Redundant nominal evaluations.** For methods whose `dx_factors` contain 0, `f(x)` is re-evaluated **once per column** (per group when partitioned). Forward differences therefore cost `2n` calls instead of `n+1` unless the cache is enabled. Consider evaluating `f(x)` once in `compute_jacobian` and reusing it. **FIXED**
4. **Bounds handling for mode 1.** When a user specifies `jacobian_method(s)`, perturbations can step outside `[xlow, xhigh]` with no warning. Only `class` mode checks bounds. That's reasonable, but worth documenting.
5. **Selection when every method violates the bounds.** `select_finite_diff_method` falls back to `meth(1)`, which is usually a central difference: the most likely method to violate a one-sided bound. The one-sided method that violates least would be a better fallback.

## Performance

- **O(n·nnz) index work.** `count(icol==i)` and `pack(..., icol==i)` scan the whole pattern for every column ([2501](https://github.com/jacobwilliams/NumDiff/blob/master/src/numerical_differentiation_module.f90#L2501), [2572](https://github.com/jacobwilliams/NumDiff/blob/master/src/numerical_differentiation_module.f90#L2572), and again inside `columns_in_partition_group`, plus the partitioned divide loop). For large sparse problems this dominates. Build a CSC pointer array (`jpntr`, which `dsm` already computes and then discards) once when the pattern is set, and slice it.
- **Quadratic growth on append.** `nonzero_rows = [nonzero_rows, …]` in `columns_in_partition_group` reallocates on every append. Precompute the group sizes, or cache the per-group index lists when the partition is built: they never change between Jacobian calls.
- `get_all_methods_in_class` is called on every `compute_sparsity_random_2` call and rebuilds the list by trying IDs one by one. Cache it once.

---

## Testing & CI

- **No assertions.** `test1`/`test2` print values and function counts but never compare against analytic derivatives, so every bug above passes CI. Suggested minimum:
  - Compare against analytic Jacobians for each method/class, dense, partitioned, and `diff`.
  - Use a tolerance based on method order.
  - Add regression tests for items 1–7 (the reproducers used for this review are short and could be adapted directly).
- **Coverage gaps.** No test covers:
  - classes 11–17
  - `sparsity_mode=2` or `3`
  - `set_sparsity_pattern` with linear elements
  - cache partial hits
  - `terminate()`
  - re-initialization
  - `REAL32`/`REAL128` builds
- **CI** (`.github/workflows/CI.yml`):
  - It pins **gfortran 10** (2020) and uses deprecated `actions/checkout@v3`/`setup-python@v4`.
  - The compile step is commented out.
  - Consider a matrix over gfortran 12–15 plus ifx.
  - Add a debug job with `-fcheck=all -ffpe-trap=invalid,zero,overflow -finit-real=snan`. That flag set alone catches #3, #5 and #8.
- **Version drift.** `fpm.toml` says `1.5.1`, while the new (untracked) `pixi.toml` says `1.8.2`.

---

## Suggested priority

1. Fix #1–#4 (wrong results or crashes in normal use) and add regression tests for them.
2. Fix #6 and #7 (a one-line fix and a two-line fix).
3. Make the hash portable and fix `divide_interval`.
4. Build the CSC index / per-group cache for performance.
5. Tidy the remaining validation and `intent` items.

# Changes to the MPL-2.0 SLIM files

Each change stays under the MPL-2.0. Dates are those of the commit that made the change. The files' copyright and licence
lines are unchanged.

- 2026-11-11: the eight files vendored unchanged (their SHA-256 in README.md).
- 2026-11-11: `slim.h`, `slim.cpp`: `slim_precompute_rest` starts the solver from a rest shape given as gradient operators
  `Dx`, `Dy` and face areas instead of vertex positions (no constraints, `V` left empty; the face weights are set to one,
  unweighted). The part of `pre_calc` after the gradients moved into `allocate_solver_state`, shared by both.
- 2026-11-11: `slim.h`, `slim.cpp`: `SLIMData::failed_iteration` (-1 while none). `solve_weighted_arap` checks the
  factorisation's and the solve's `info()` and that every value is finite; on a failure it sets `failed_iteration` and
  leaves the map unchanged instead of throwing, and `slim_solve` stops its loop, so `V_o` stays the last map the
  flip-avoiding line search accepted.
- 2026-11-11: `flip_avoiding_line_search.cpp`: `get_smallest_pos_quad_zero` with `a == 0` (the determinant `b t + c` linear in
  the step) returns the root `t = -c / b` when it is positive and finite, and infinity otherwise or when `b == 0`. It took
  `-b / c`, so it never limited the step in that case and threw `SlimFailedException` when `c == 0`.
- 2026-11-12: all eight files: the Doxygen comment block below the licence header (`\file` and `\ingroup`) is removed.
- 2026-11-12: `slim.cpp`, `doublearea.cpp`, `edge_lengths.cpp`: the three stand-in headers have new names, so the include
  lines read `slim_assert.h`, `slim_math_base.h` and `slim_parallel.hh`; the stand-ins' namespace is `slim_support`, so
  `using namespace` names it; the assert macro is `SLIM_ASSERT`. Nothing else changes: the same arithmetic in the same
  order, and the assert checks nothing as before.

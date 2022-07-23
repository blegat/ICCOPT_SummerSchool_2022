# ICCOPT 2022 JuMP Tutorial

This repository will contain materials for the Julia/JuMP tutorial at
[ICCOPT 2022](https://iccopt2022.lehigh.edu/summer-school/summer-school-program/).

## Installation requirements

Participants are expected to have their own laptop to follow along with the tutorial.

Please install:
- Julia v1.6 or v1.7 from https://julialang.org/downloads/.
- [Visual Studio Code (VSCode)](https://code.visualstudio.com/) (text editor) with [Julia extensions](https://code.visualstudio.com/docs/languages/julia).
- The `jupyter` VSCode extension and IJulia with `julia -e 'using Pkg; Pkg.add("IJulia"); Pkg.build("IJulia");'`.

### Common issues

See Julia's
[platform-specific instructions](https://julialang.org/downloads/platform/#platform_specific_instructions_for_official_binaries)
for more guidance.

If you have trouble installing before the tutorial, please
[open a GitHub issue](https://github.com/blegat/ICCOPT_SummerSchool_2022/issues/new) and we'll try to help.

## Download these materials

**Thanks for being prepared and checking out the materials before the tutorial. We're still working on them; please wait for this message to disappear before downloading.**

Next, download a copy of these materials.

**If you have `git` installed**

`cd` to an appropriate directory, then run
```
git clone https://github.com/blegat/ICCOPT_SummerSchool_2022
```

**If you don't have `git` installed**

[Download this zip file](https://github.com/blegat/ICCOPT_SummerSchool_2022/archive/main.zip).
Once downloaded, unzip it to an appropriate location.

## Extra resources

- [Julia's documentation](https://docs.julialang.org/en/v1/)
- [JuMP's documentation](https://jump.dev/JuMP.jl/stable/)
- [Discourse forum](https://discourse.julialang.org/c/domain/opt/13)
- [JuliaCon 2022](https://juliacon.org/2022/) and [JuMP-dev 2022](https://jump.dev/meetings/juliacon2022/)

## Schedule

### Environment setup and basics of the Julia programming language

Time: 09:00 – 10:30

Open the notebook `getting_started_with_julia.ipynb`.

### Introduction to JuMP

Time: 11:00 – 12:30

Open the notebook `getting_started_with_JuMP.ipynb`.

### Modeling nonlinear and conic optimization problems with JuMP

Time: 02:00 – 03:30

Open the notebooks:
* `passport_problem.ipynb` then
* `portfolio_optimization.ipynb` then
* `rocket_control.ipynb`.

### How to write your own conic solver

Time: 04:00 – 05:30

Open the REPL with `Ctrl+Shift+P` (in VSCode) then write `Julia: Start REPL`:
```julia
julia> using Revise

julia> pwd()
"/.../ICCOPT_SummerSchool_2022"

julia> using Pkg

julia> Pkg.activate(".")
  Activating project at `/.../ICCOPT_SummerSchool_2022`

julia> Pkg.develop(path="SimpleConicADMM")

julia> include("SimpleConicADMM/test/JuMP.jl")

Iter. | Primal Feas.   | Primal Obj.    | Dual Feas.     | Dual Obj.      | Rel. gap    | Time (s)
    0 |  +5.000000e-01 |  +0.000000e+00 |  +5.857864e-01 |  -0.000000e+00 |   0.000e+00 |   5.364e-01
    1 |  +5.000000e-01 |  +0.000000e+00 |  +5.857864e-01 |  -0.000000e+00 |   0.000e+00 |   1.518e+00
(...)
```
Now, open the file `SimpleConicADMM/src/solver.jl` and fix the `FIXME`s and do
the `TODO`s (follow equation (17) of the
[SCS paper](https://web.stanford.edu/~boyd/papers/pdf/scs_long.pdf)) and then try again:
```julia
julia> test_SOC()

Iter. | Primal Feas.   | Primal Obj.    | Dual Feas.     | Dual Obj.      | Rel. gap    | Time (s)
    0 |  +5.000000e-01 |  +0.000000e+00 |  +5.857864e-01 |  -0.000000e+00 |   0.000e+00 |   5.364e-01
    1 |  +5.000000e-01 |  +0.000000e+00 |  +5.857864e-01 |  -0.000000e+00 |   0.000e+00 |   1.518e+00
(...)
```
It should take your changes into account thanks to `Revise` unless you change
the definition of a `struct`, e.g., add a field to `struct Cache`.
In that case close the REPL (`Ctrl+D`), and start again.

# Lifetime Exposures and Asthma outcomes Projection (LEAP)

Documentation for LEAP will be available by the end of 2023.

Here is the pseudo-code of the model:

![Pseudocode](documentation/pseudocode.png)

# Developers

## Installation

To install this package on your computer locally, first download it from `GitHub`:

```
git clone https://github.com/resplab/LEAP.jl
```

Next, install the package:

```
cd LEAP.jl
julia
> import Pkg
> Pkg.add(path="PATH_TO/LEAP.jl")
```

## Dependencies

To add dependencies, you can add them to the `LEAP.jl` file:

```
using MyDependency
```

For the dependencies to be included in the package, you need to also add them to the `Project.toml`
file. To do so, open a terminal and enter the `Julia` package `REPL`:

```
cd LEAP.jl
julia
julia> # type ] here to enter the Pkg REPL
[(@v1.x) pkg> activate . # activate our package REPL
[(LEAP) pkg> add MyDependency
```

## Tests

### Command Line

To run unit tests via the command line:

```
cd LEAP.jl
julia --project=. --code-coverage test/runtests.jl
```

To select specific modules to test, list them as args:

```
julia --project=. --code-coverage test/runtests.jl "agent"
```

If you have run with `--code-coverage`, then after the tests have run, you will need to run
the `teardown.jl` file:

```
julia --project=. "./test/teardown.jl" --module-coverage "agent"
```


### Julia REPL

```
cd LEAP.jl
julia
julia> using Pkg
julia> Pkg.add(path="PATH_TO/LEAP.jl")
julia> Pkg.activate(@__DIR__)
julia> Pkg.instantiate()
julia> Pkg.test("LEAP", coverage=true)
```

To select specific modules to test, use the `test_args` argument:

```
julia> Pkg.test("LEAP", coverage=true, tests_args=["all"])
```

### Package Mode

```
cd LEAP.jl
julia
julia> # type ] here to enter the Pkg REPL
[(@v1.x) pkg> activate . # activate our package REPL
[(LEAP) pkg> test
```

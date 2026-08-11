# API Reference

```@index
```

```@autodocs
Modules = [ 
        H2Trees,
        # `@autodocs` does not recurse into submodules, so the internal Hilbert submodule is
        # listed explicitly; otherwise `checkdocs` reports its docstrings as undocumented.
        H2Trees.HilbertOrdering,
        H2Trees.SEBB,
        if isdefined(Base, :get_extension)
            Base.get_extension(H2Trees, :H2MetisTrees)
        else
            H2Trees.H2MetisTrees
        end,
        if isdefined(Base, :get_extension)
            Base.get_extension(H2Trees, :H2BEASTTrees)
        else
            H2Trees.H2BEASTTrees
        end,
        if isdefined(Base, :get_extension)
            Base.get_extension(H2Trees, :H2ParallelKMeansTrees)
        else
            H2Trees.H2ParallelKMeansTrees
        end,
        if isdefined(Base, :get_extension)
            Base.get_extension(H2Trees, :H2PlotlyJSTrees)
        else
            H2Trees.H2PlotlyJSTrees
        end,
        ]
```

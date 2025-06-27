# API Reference

```@index
```

```@autodocs
Modules = [ 
        H2Trees,
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

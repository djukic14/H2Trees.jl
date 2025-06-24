using Test
using BEAST, CompScienceMeshes
using StaticArrays
using LinearAlgebra
using H2Trees

@testset "AllTranslations TwoNTree" begin
    λ = 1.0
    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere6.in")
    )
    X = raviartthomas(m)
    tree = TwoNTree(X, λ / 10; minvalues=20)

    disaggregationplan = H2Trees.DisaggregateTranslatePlan(
        tree, H2Trees.TranslatingNodesIterator
    )
    aggregationplan = H2Trees.AggregateTranslatePlan(tree, H2Trees.TranslatingNodesIterator)

    translationinfo, translations = H2Trees.translations(
        tree, disaggregationplan, H2Trees.AllTranslations()
    )
    atranslationinfo, atranslations = H2Trees.translations(
        tree, aggregationplan, H2Trees.AllTranslations()
    )

    numberoftranslationinfo = 0
    for i in eachindex(translationinfo)
        numberoftranslationinfo += length(translationinfo[i])
    end
    @test numberoftranslationinfo == length(translations)

    tfs = Dict{NTuple{2,Int},eltype(translations)}()

    for i in eachindex(translationinfo)
        for j in eachindex(translationinfo[i])
            tfs[(translationinfo[i][j].receivingnode, translationinfo[i][j].translatingnode)] = translations[translationinfo[i][j].translationID]
        end
    end

    atfs = Dict{NTuple{2,Int},eltype(translations)}()
    for i in eachindex(atranslationinfo)
        for j in eachindex(atranslationinfo[i])
            atfs[(atranslationinfo[i][j].receivingnode, atranslationinfo[i][j].translatingnode)] = atranslations[atranslationinfo[i][j].translationID]
        end
    end

    for receivingnode in H2Trees.DepthFirstIterator(tree, 1)
        for translatingnode in H2Trees.TranslatingNodesIterator(tree, receivingnode)
            @test tfs[(receivingnode, translatingnode)] ≈
                H2Trees.center(tree, receivingnode) -
                  H2Trees.center(tree, translatingnode)
            @test atfs[(receivingnode, translatingnode)] ≈
                tfs[(receivingnode, translatingnode)]
        end
    end
end

@testset "DirectionInvariance(PerLevel) TwoNTree" begin
    λ = 1.0
    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere6.in")
    )
    X = raviartthomas(m)
    tree = TwoNTree(X, λ / 10; minvalues=30)

    disaggregationplan = H2Trees.DisaggregateTranslatePlan(
        tree, H2Trees.TranslatingNodesIterator
    )
    aggregationplan = H2Trees.AggregateTranslatePlan(tree, H2Trees.TranslatingNodesIterator)

    for TranslationTrait in
        [H2Trees.DirectionInvariancePerLevel(), H2Trees.DirectionInvariance()]
        translationinfo, translations = H2Trees.translations(
            tree, disaggregationplan, TranslationTrait
        )
        atranslationinfo, atranslations = H2Trees.translations(
            tree, aggregationplan, TranslationTrait
        )

        numberoftranslationinfo = 0
        for i in eachindex(translationinfo)
            numberoftranslationinfo += length(translationinfo[i])
        end
        @test numberoftranslationinfo != length(translations)

        tfs = Dict{NTuple{2,Int},eltype(translations)}()
        for i in eachindex(translationinfo)
            for j in eachindex(translationinfo[i])
                tfs[(translationinfo[i][j].receivingnode, translationinfo[i][j].translatingnode)] = translations[translationinfo[i][j].translationID]
            end
        end

        atfs = Dict{NTuple{2,Int},eltype(translations)}()
        for i in eachindex(atranslationinfo)
            for j in eachindex(atranslationinfo[i])
                atfs[(atranslationinfo[i][j].receivingnode, atranslationinfo[i][j].translatingnode)] = atranslations[atranslationinfo[i][j].translationID]
            end
        end

        for receivingnode in H2Trees.DepthFirstIterator(tree, 1)
            for translatingnode in H2Trees.TranslatingNodesIterator(tree, receivingnode)
                @test tfs[(receivingnode, translatingnode)] ≈
                    H2Trees.center(tree, receivingnode) -
                      H2Trees.center(tree, translatingnode)
                @test atfs[(receivingnode, translatingnode)] ≈
                    tfs[(receivingnode, translatingnode)]
            end
        end
    end
end

@testset "AllTranslations Block-TwoNTree" begin
    λ = 1.0
    mx = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere6.in")
    )
    my = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "cuboid2.in")
    )
    X = raviartthomas(mx)
    Y = raviartthomas(my)
    tree = TwoNTree(Y, X, λ / 10; minvaluestest=20, minvaluestrial=20)

    disaggregationplan = H2Trees.DisaggregateTranslatePlan(
        H2Trees.testtree(tree), H2Trees.trialtree(tree), H2Trees.TranslatingNodesIterator
    )
    aggregationplan = H2Trees.AggregateTranslatePlan(
        H2Trees.trialtree(tree), H2Trees.testtree(tree), H2Trees.TranslatingNodesIterator
    )

    translationinfo, translations = H2Trees.translations(
        tree, disaggregationplan, H2Trees.AllTranslations()
    )
    atranslationinfo, atranslations = H2Trees.translations(
        tree, aggregationplan, H2Trees.AllTranslations()
    )

    numberoftranslationinfo = 0
    for i in eachindex(translationinfo)
        numberoftranslationinfo += length(translationinfo[i])
    end
    @test numberoftranslationinfo == length(translations)

    tfs = Dict{NTuple{2,Int},eltype(translations)}()

    for i in eachindex(translationinfo)
        for j in eachindex(translationinfo[i])
            tfs[(translationinfo[i][j].receivingnode, translationinfo[i][j].translatingnode)] = translations[translationinfo[i][j].translationID]
        end
    end

    atfs = Dict{NTuple{2,Int},eltype(translations)}()
    for i in eachindex(atranslationinfo)
        for j in eachindex(atranslationinfo[i])
            atfs[(atranslationinfo[i][j].receivingnode, atranslationinfo[i][j].translatingnode)] = atranslations[atranslationinfo[i][j].translationID]
        end
    end

    for receivingnode in H2Trees.DepthFirstIterator(H2Trees.testtree(tree), 1)
        for translatingnode in H2Trees.TranslatingNodesIterator(
            H2Trees.trialtree(tree), H2Trees.testtree(tree), receivingnode
        )
            @test tfs[(receivingnode, translatingnode)] ==
                H2Trees.center(H2Trees.testtree(tree), receivingnode) -
                  H2Trees.center(H2Trees.trialtree(tree), translatingnode)
            @test atfs[(receivingnode, translatingnode)] ≈
                tfs[(receivingnode, translatingnode)]
        end
    end
end

@testset "DirectionInvariance(PerLevel) Block-TwoNTree" begin
    λ = 1.0
    mx = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere6.in")
    )
    my = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "cuboid2.in")
    )

    X = raviartthomas(mx)
    Y = raviartthomas(my)
    tree = TwoNTree(X, Y, λ / 10; minvaluestest=20, minvaluestrial=20)

    disaggregationplan = H2Trees.DisaggregateTranslatePlan(
        H2Trees.testtree(tree), H2Trees.trialtree(tree), H2Trees.TranslatingNodesIterator
    )
    aggregationplan = H2Trees.AggregateTranslatePlan(
        H2Trees.trialtree(tree), H2Trees.testtree(tree), H2Trees.TranslatingNodesIterator
    )

    for TranslationTrait in
        [H2Trees.DirectionInvariancePerLevel(), H2Trees.DirectionInvariance()]
        translationinfo, translations = H2Trees.translations(
            tree, disaggregationplan, TranslationTrait
        )
        atranslationinfo, atranslations = H2Trees.translations(
            tree, aggregationplan, TranslationTrait
        )

        numberoftranslationinfo = 0
        for i in eachindex(translationinfo)
            numberoftranslationinfo += length(translationinfo[i])
        end
        @test numberoftranslationinfo != length(translations)

        tfs = Dict{NTuple{2,Int},eltype(translations)}()
        for i in eachindex(translationinfo)
            for j in eachindex(translationinfo[i])
                tfs[(translationinfo[i][j].receivingnode, translationinfo[i][j].translatingnode)] = translations[translationinfo[i][j].translationID]
            end
        end

        atfs = Dict{NTuple{2,Int},eltype(translations)}()
        for i in eachindex(atranslationinfo)
            for j in eachindex(atranslationinfo[i])
                atfs[(atranslationinfo[i][j].receivingnode, atranslationinfo[i][j].translatingnode)] = atranslations[atranslationinfo[i][j].translationID]
            end
        end

        for receivingnode in H2Trees.DepthFirstIterator(H2Trees.testtree(tree), 1)
            for translatingnode in H2Trees.TranslatingNodesIterator(
                H2Trees.trialtree(tree), H2Trees.testtree(tree), receivingnode
            )
                @test tfs[(receivingnode, translatingnode)] ≈
                    H2Trees.center(H2Trees.testtree(tree), receivingnode) -
                      H2Trees.center(H2Trees.trialtree(tree), translatingnode)
                @test atfs[(receivingnode, translatingnode)] ≈
                    tfs[(receivingnode, translatingnode)]
            end
        end
    end
end

@testset "DirectionInvariancePerLevel Comparison" begin
    # there are two version for DirectionInvariancePerLevel: one general and one for TwoNTree
    # here we compare both methods

    λ = 1.0
    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere6.in")
    )
    X = raviartthomas(m)
    tree = TwoNTree(X, λ / 10; minvalues=30)

    disaggregationplan = H2Trees.DisaggregateTranslatePlan(
        tree, H2Trees.TranslatingNodesIterator
    )
    aggregationplan = H2Trees.AggregateTranslatePlan(tree, H2Trees.TranslatingNodesIterator)

    TranslationTraits = [
        H2Trees.DirectionInvariancePerLevel(),
        H2Trees.DirectionInvariance(),
        H2Trees.AllTranslations(),
    ]
    for TranslationTrait in TranslationTraits
        translationinfo, translations = H2Trees.translations(
            tree, disaggregationplan, TranslationTrait
        )
        translationinfogeneral, translationsgeneral = H2Trees.translations(
            tree, H2Trees.isAnyTree(), disaggregationplan, TranslationTrait
        )
        atranslationinfogeneral, atranslationsgeneral = H2Trees.translations(
            tree, H2Trees.isAnyTree(), aggregationplan, TranslationTrait
        )

        tfs = Dict{NTuple{2,Int},eltype(translations)}()
        for i in eachindex(translationinfo)
            for j in eachindex(translationinfo[i])
                tfs[(translationinfo[i][j].receivingnode, translationinfo[i][j].translatingnode)] = translations[translationinfo[i][j].translationID]
            end
        end

        atfs = Dict{NTuple{2,Int},eltype(translations)}()
        for i in eachindex(atranslationinfogeneral)
            for j in eachindex(atranslationinfogeneral[i])
                atfs[(atranslationinfogeneral[i][j].receivingnode, atranslationinfogeneral[i][j].translatingnode)] = atranslationsgeneral[atranslationinfogeneral[i][j].translationID]
            end
        end

        tfsgeneral = Dict{NTuple{2,Int},eltype(translationsgeneral)}()
        for i in eachindex(translationinfogeneral)
            for j in eachindex(translationinfogeneral[i])
                tfsgeneral[(translationinfogeneral[i][j].receivingnode, translationinfogeneral[i][j].translatingnode)] = translationsgeneral[translationinfogeneral[i][j].translationID]
            end
        end

        for receivingnode in H2Trees.DepthFirstIterator(tree, 1)
            for translatingnode in H2Trees.TranslatingNodesIterator(tree, receivingnode)
                @test tfs[(receivingnode, translatingnode)] ≈
                    tfsgeneral[(receivingnode, translatingnode)] ≈
                    H2Trees.center(tree, receivingnode) -
                    H2Trees.center(tree, translatingnode) ≈
                    atfs[(receivingnode, translatingnode)]
            end
        end
    end

    λ = 1.0
    mx = m
    X = raviartthomas(mx)

    for my in [
        CompScienceMeshes.readmesh(
            joinpath(pkgdir(H2Trees), "test", "assets", "in", "spherewithcenter11.in")
        ),
        CompScienceMeshes.readmesh(
            joinpath(pkgdir(H2Trees), "test", "assets", "in", "spherewithcenter12.in")
        ),
        CompScienceMeshes.readmesh(
            joinpath(pkgdir(H2Trees), "test", "assets", "in", "spherewithcenter13.in")
        ),
    ]
        for TranslationTrait in TranslationTraits
            Y = raviartthomas(my)
            tree2 = TwoNTree(X, Y, λ / 10; minvaluestest=20, minvaluestrial=30)

            disaggregationplan2 = H2Trees.DisaggregateTranslatePlan(
                H2Trees.testtree(tree2),
                H2Trees.trialtree(tree2),
                H2Trees.TranslatingNodesIterator,
            )

            aggregationplan2 = H2Trees.AggregateTranslatePlan(
                H2Trees.trialtree(tree2),
                H2Trees.testtree(tree2),
                H2Trees.TranslatingNodesIterator,
            )

            translationinfo2, translations2 = H2Trees.translations(
                tree2, disaggregationplan2, TranslationTrait
            )

            translationinfogeneral2, translationsgeneral2 = H2Trees.translations(
                H2Trees.testtree(tree2),
                H2Trees.trialtree(tree2),
                H2Trees.isAnyTree(),
                H2Trees.isAnyTree(),
                disaggregationplan2,
                TranslationTrait,
            )

            atranslationinfogeneral2, atranslationsgeneral2 = H2Trees.translations(
                H2Trees.testtree(tree2),
                H2Trees.trialtree(tree2),
                H2Trees.isAnyTree(),
                H2Trees.isAnyTree(),
                aggregationplan2,
                TranslationTrait,
            )

            tfs2 = Dict{NTuple{2,Int},eltype(translations2)}()
            for i in eachindex(translationinfo2)
                for j in eachindex(translationinfo2[i])
                    tfs2[(translationinfo2[i][j].receivingnode, translationinfo2[i][j].translatingnode)] = translations2[translationinfo2[i][j].translationID]
                end
            end

            tfsgeneral2 = Dict{NTuple{2,Int},eltype(translationsgeneral2)}()
            for i in eachindex(translationinfogeneral2)
                for j in eachindex(translationinfogeneral2[i])
                    tfsgeneral2[(translationinfogeneral2[i][j].receivingnode, translationinfogeneral2[i][j].translatingnode)] = translationsgeneral2[translationinfogeneral2[i][j].translationID]
                end
            end

            atfsgeneral2 = Dict{NTuple{2,Int},eltype(translationsgeneral2)}()
            for i in eachindex(atranslationinfogeneral2)
                for j in eachindex(atranslationinfogeneral2[i])
                    atfsgeneral2[(atranslationinfogeneral2[i][j].receivingnode, atranslationinfogeneral2[i][j].translatingnode)] = atranslationsgeneral2[atranslationinfogeneral2[i][j].translationID]
                end
            end

            for receivingnode in H2Trees.DepthFirstIterator(H2Trees.testtree(tree2), 1)
                for translatingnode in H2Trees.TranslatingNodesIterator(
                    H2Trees.trialtree(tree2), H2Trees.testtree(tree2), receivingnode
                )
                    @test tfs2[(receivingnode, translatingnode)] ≈
                        tfsgeneral2[(receivingnode, translatingnode)] ≈
                        H2Trees.center(H2Trees.testtree(tree2), receivingnode) -
                        H2Trees.center(H2Trees.trialtree(tree2), translatingnode) ≈
                        atfsgeneral2[(receivingnode, translatingnode)]
                end
            end
        end
    end
end

@testset "AllLeavesTranslationsIterator" begin
    λ = 1.0
    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere6.in")
    )

    tree = TwoNTree(raviartthomas(m), λ / 10; minvalues=10)

    f = H2Trees.AllLeavesTranslationsIterator(tree)

    for leaf in H2Trees.leaves(tree)
        level = H2Trees.level(tree, leaf)
        comparison = sort(collect(H2Trees.leaves(tree)))
        comparison = filter(x -> H2Trees.level(tree, x) == level, comparison)
        @test sort(collect(f(tree, leaf))) == sort(comparison)
    end
end

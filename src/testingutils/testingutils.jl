
function testwellseparatedness(tree)
    return testwellseparatedness(tree, treetrait(tree))
end

function testwellseparatedness(tree, ::Any)
    @threads for level in levels(tree)
        for trialnode in LevelIterator(tree, level)
            for testnode in LevelIterator(tree, level)
                testnode == trialnode && continue

                connections = 0

                neartranslations = 0
                directtranslations = 0
                parenttranslations = 0

                if isnear(tree, testnode, trialnode, treetrait(tree))
                    connections += 1
                    neartranslations += 1
                end

                if testnode in TranslatingNodesIterator(tree, trialnode)
                    connections += 1
                    directtranslations += 1
                end

                trialparent = trialnode

                for testparent in ParentUpwardsIterator(tree, testnode)
                    testparent == 0 && break

                    trialparent = parent(tree, trialparent)

                    if testparent in TranslatingNodesIterator(tree, trialparent)
                        connections += 1
                        parenttranslations += 1
                    end
                end

                if connections != 1
                    println(
                        "Detected connection error: parenttranslations: $parenttranslations directtranslations: $directtranslations neartranslations: $neartranslations for testnode $testnode and trialnode $trialnode on level $level",
                    )

                    return false
                end
            end
        end
    end
    return true
end

function testwellseparatedness(tree, ::isBlockTree)
    trialtree = H2Trees.trialtree(tree)
    testtree = H2Trees.testtree(tree)

    for level in levels(trialtree)
        for trialnode in LevelIterator(trialtree, level)
            for testnode in LevelIterator(testtree, level)
                connections = 0
                directranslations = 0
                parentranslations = 0

                if isnear(testtree, trialtree, testnode, trialnode)
                    connections += 1
                    directranslations += 1
                end

                if testnode in TranslatingNodesIterator(testtree, trialtree, trialnode)
                    connections += 1
                    directranslations += 1
                end

                trialparent = trialnode

                for testparent in ParentUpwardsIterator(testtree, testnode)
                    testparent == 0 && break

                    trialparent = parent(trialtree, trialparent)
                    trialparent == 0 && break

                    if testparent in
                        TranslatingNodesIterator(testtree, trialtree, trialparent)
                        connections += 1
                        parentranslations += 1
                    end
                end
                if connections != 1
                    println(
                        "Detected connection error: parenttranslations: $parenttranslations directtranslations: $directtranslations neartranslations: $neartranslations for testnode $testnode and trialnode $trialnode on level $level",
                    )

                    return false
                end
            end
        end
    end
    return true
end

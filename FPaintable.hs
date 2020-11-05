module FPaintable where

import Control.Monad(filterM)

type Parts = [Int]

type FVals = [[Int]]

check :: Parts -> FVals -> Bool
check ps fs = all (\(p,f) -> p == length f) $ zip ps fs

subsets :: Int -> [[Int]]
subsets n = filterM (\x -> [True, False]) [1..n]

lister :: Parts -> [[[Int]]]
lister ps = map subsets ps -- not quite right

fPaintable :: Parts -> FVals -> Maybe Bool
fPaintable ps fs = case check ps fs of
                        False -> Nothing
                        True -> Just True

module FPaintable where

import Control.Monad(filterM)

type Parts = [Int]

type Choice = [[Int]]

type FVals = [[Int]]

check :: Parts -> FVals -> Bool
check ps fs = all (\(p,f) -> p == length f) $ zip ps fs

subsets :: Int -> [[Int]]
subsets n = filterM (\x -> [False, True]) [1..n]

lister :: Parts -> [Choice]
lister ps = filter (any (/=[])) $ sequence $ map subsets ps

painter :: Parts -> FVals -> Choice -> [(Parts, FVals)]
painter ps fs cs = []

fPaintable' :: Parts -> FVals -> Bool
fPaintable' [] [] = True
fPaintable' (0:ps) ([]:fs) = fPaintable' ps fs
fPaintable' ps fs = all (any (uncurry fPaintable') . painter ps fs) $ lister ps

fPaintable :: Parts -> FVals -> Maybe Bool
fPaintable [] [] = Just True
fPaintable (0:ps) ([]:fs) = fPaintable ps fs
fPaintable ps fs = case check ps fs of
                        False -> Nothing
                        True -> Just $ all (any (uncurry fPaintable') . painter ps fs) $ lister ps

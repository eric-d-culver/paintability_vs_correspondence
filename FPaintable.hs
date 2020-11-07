module FPaintable(fPaintable, paintability) where

import Control.Monad(filterM)

type Parts = [Int]

type Choice = [[Int]]

type FVals = [[Int]]

check :: Parts -> FVals -> Bool
check ps fs = all (\(p,f) -> p == length f) $ zip ps fs

subsets :: Int -> [[Int]]
subsets n = filterM (\x -> [False, True]) [0..(n-1)]

lister :: Parts -> [Choice]
lister ps = filter (any (/=[])) $ sequence $ map subsets ps

new_parts :: Parts -> Choice -> Int -> Parts
new_parts ps cs n = let (ys, zs) = splitAt n ps in ys ++ [(head zs) - (length (cs !! n))] ++ (tail zs)
--new_parts (p:ps) (c:cs) 0 = (p - (length c)) : ps
--new_parts (p:ps) (c:cs) n = p : (new_parts ps cs (n-1))
--new_parts [] _ _ = []

delete_index :: Int -> [a] -> [a]
delete_index n xs = let (ys, zs) = splitAt n xs in ys ++ (tail zs)

dec_index :: (Num a) => Int -> [a] -> [a]
dec_index n xs = let (ys, zs) = splitAt n xs in ys ++ [(head zs) - 1] ++ (tail zs)

compose :: [a -> a] -> a -> a
compose fs v = foldl (flip (.)) id fs $ v

new_fvals :: FVals -> Choice -> Int -> FVals
new_fvals (f:fs) (c:cs) (-1) = (compose (map dec_index c) f) : (new_fvals fs cs (-1))
new_fvals (f:fs) (c:cs) 0 = (compose (map delete_index c) f) : (new_fvals fs cs (-1))
new_fvals (f:fs) (c:cs) n = (compose (map dec_index c) f) : (new_fvals fs cs (n-1))
new_fvals [] [] _ = []

new_pairs :: Parts -> FVals -> Choice -> Int -> (Parts, FVals)
new_pairs ps fs cs i = (new_parts ps cs i, new_fvals fs cs i)

painter :: Parts -> FVals -> Choice -> [(Parts, FVals)]
painter ps fs cs = map (new_pairs ps fs cs) $ filter (\i -> cs !! i /= []) [0..(length ps - 1)]

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

all_k :: Parts -> Int -> FVals
all_k ps k = map (flip replicate k) ps

paintability :: Parts -> Int
paintability ps = head $ filter (fPaintable' ps . all_k ps) [2..]

module FPaintable(fPaintable, paintability) where

import Control.Monad(filterM)

type Parts = [Int]

type Choice = [[Int]]

type FVals = [[Int]]

leqList :: [Int] -> [[Int]]
leqList l = sequence $ map (\i -> [0..i]) l

lister :: FVals -> [Choice]
lister fs = filter (any (\x -> sum(x) /= 0)) $ sequence $ map leqList fs 

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
-- could sort each list, might make it faster because of symmetries

painter :: FVals -> Choice -> [FVals]
painter fs cs = map (new_fvals fs cs) $ filter (\i -> (cs !! i) /= []) [0..(length fs - 1)]

fPaintable :: FVals -> Bool
fPaintable [] = True
fPaintable ([]:fs) = fPaintable fs
fPaintable fs = all (any fPaintable . painter fs) $ lister fs

all_k :: Parts -> Int -> FVals
all_k ps k = map (flip replicate k) ps

paintability :: Parts -> Int
paintability ps = head $ filter (fPaintable . all_k ps) [2..]

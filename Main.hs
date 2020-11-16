import FPaintable

main = do 
        putStrLn $ show $ fPaintable [2,2] [[2,2],[2,2]]
        putStrLn $ show $ fPaintable [4,4] [[2,2,2,2],[2,2,2,2]]
        putStrLn $ show $ fPaintable [4,4] [[3,3,3,3],[3,3,3,3]]
        putStrLn $ show $ paintability [3,3]
        putStrLn $ show $ paintability [4,4]

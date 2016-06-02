import Statistics.Sample.KernelDensity
import Data.Vector.Unboxed as V
import Graphics.Gnuplot.Simple
import System.Environment
import Data.List
    
process :: Vector Double -> Int -> IO ()
process v n = let e = kde n v
                  pts = toList $ V.zip (fst e) (snd e) 
              in  sequence_ $ fmap (\x -> putStr ((Prelude.++) (show x) " ")) pts

--plotIt :: Vector Double -> Int -> IO ()
--plotIt v n = let e = kde n v
--                 pts = toList $ V.zip (fst e) (snd e)

                -- plotList [] pts
                      
readF file n = do
  w <- readFile file
  process (fromList (stringToList w)) n

stringToList :: String -> [Double]
stringToList s =
    Data.List.reverse $ Data.List.drop 1 $ Data.List.reverse $ (Data.List.drop 1) $ fmap (\x -> read x :: Double) $ fmap (Prelude.drop 1) $ groupBy (\x y -> y /= '\n') ((Prelude.++) "\n" s)
        
        
main :: IO ()
main = do
  [file,n] <- getArgs
  readF file (read n :: Int)


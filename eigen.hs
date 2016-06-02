import Numeric.LinearAlgebra
import Data.Ratio
import Control.Monad as M
import Data.Complex
import System.Environment
    
coprime a b = (gcd a b == 1)
    
-- Generate all Farey fractions with denominator = q

farey :: Integer -> [Rational]
farey q = fmap (\x -> (x % q)) $ Prelude.filter (coprime q) [1..q]

-- Generate all Farey fractions with denominator < q

allFarey :: Integer -> [Rational]
allFarey q = M.join $ fmap (farey) [1..q]
             
sizeFarey = length . allFarey
             
-- Generate the matrix A of dimension q n (with q Farey Fractions and n the integers)

genSquareMatrix :: Integer -> Int -> Matrix (Complex Double)
genSquareMatrix q n =
    let d = sizeFarey q
    in 
      (n><d) $ M.join $ fmap (\x -> fmap (trans x) (allFarey q)) [1..n]
          where
            trans x = (mkPolar 1) . ((*) (2*pi)) . fromRational . ((*) (toRational x))
        
-- Generate the actual matrix

eigenLargeSieve q n = let d = fromIntegral (sizeFarey q)
                      in  fmap (\x -> (x*x) / d) $ toList $ singularValues (genSquareMatrix q n)

main :: IO ()
main = do
  [q,n] <- getArgs
  putStrLn ("q = " ++ q ++ ", n = " ++ n)
  sequence_ $ reverse $ fmap (putStrLn . show) (eigenLargeSieve (read q :: Integer) (read n :: Int))
                               

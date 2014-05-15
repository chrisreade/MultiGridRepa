{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE TemplateHaskell     #-}
{-# LANGUAGE QuasiQuotes         #-}
{-# LANGUAGE ScopedTypeVariables #-}

{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}


import Data.Array.Repa				    as R

import MultiGrid_Poisson

import Data.Array.Repa.IO.Timing
import System.Environment
import Data.Array.Repa.Algorithms.Pixel
import Data.Array.Repa.Algorithms.ColorRamp
import Data.Array.Repa.IO.BMP	


main :: IO ()
main 
 = do	args	<- getArgs
	case args of
	  [stackSize, filename]	
	    -> 	poissonTest (read stackSize) filename

	  _ -> do
		putStr usage
		return ()

-- | Command line usage information.
usage	:: String
usage	= unlines
	[ "Usage: poissonTest <stackSize> <outFileName>"
	, ""
	, "  stackSize  :: Int  Number of Grids to use."
	, "  outFileName  :: String  (BMP File name)"
	, "" ]

poissonTest :: Int -> String -> IO()
poissonTest gridStackSize filename =
 let 
     shapeInit = fineGridShape gridStackSize
     h  = hInit shapeInit
     f  = fInit shapeInit 
     bM = boundMask shapeInit
     bV = boundValues shapeInit
  in
    do 
      (arrFinal, t) <- time $ fullMG h f bM bV
      putStr (prettyTime t)
      writeHeatMapBMP filename arrFinal


writeHeatMapBMP :: String 
                -> Array U DIM2 Double
                -> IO()
writeHeatMapBMP !filename !arr
  = do !maxVal <- foldAllP max 0.0 arr
       !minVal <- foldAllP min 0.0 arr
       !arrImageOut <- computeP
                       $  R.map rgb8OfDouble
                       $  R.map (rampColorHotToCold minVal maxVal) arr
       writeImageToBMP filename arrImageOut
    

fineGridShape:: Int -> DIM2
fineGridShape gridStackSize 
  = (Z :. n :. n) where n = 2^gridStackSize+1


{- Example 2 https://math.la.asu.edu/~kuiper/502files/Laplace.pdf Exercise 11
    
    >>>> Au = -f <<<< 

   u and f defined on the square (0,0) to (pi,pi) (u is unknown)

   f(x,y) = -1

   Dirichlet Boundary Conditions

   u(x,0) = 0 and u(x,pi) = 0            for 0<=x<=pi
   u(0,y) = 0                            for 0<=y<=pi
   u(pi,y) = y*(pi-y)                    for 0<=y<=pi
-}

distance:: Double
distance = pi  --length of sides for square area covered by grids

hInit :: DIM2 -> Double
hInit sh = 
  case sh of
    (Z:.m:.n) | m==n      -> distance / fromIntegral (m-1) 
    _         | otherwise -> error "non-square grid shape"

hCoords :: DIM2 -> [Double]
hCoords sh
  = Prelude.map ((h*) . fromIntegral) [0..n-1]
    where
        (Z :. _ :. n) = sh
        h = hInit sh

boundMask :: DIM2 -> Array U DIM2 Double             
boundMask sh = 
  fromListUnboxed sh $ concat 
     [
      take ((m-1)*n) (cycle mainRow),
      edgeRow
     ]
     where (Z:.m:.n) = sh
           mainRow = 0.0: replicate (n-2) 1.0  Prelude.++ [0.0]
           edgeRow = replicate n 0.0


boundValues :: DIM2 -> Array U DIM2 Double        
boundValues sh =
  fromListUnboxed sh $ concat
     [
       take ((m-1)*n) (cycle mainRow),
       edgeRow
     ] 
     where (Z:.m:.n) = sh
           mainRow = replicate n 0.0
           edgeRow = Prelude.map (\j -> j*(pi-j)) coords
           coords = hCoords sh

fInit :: DIM2 -> Array U DIM2 Double
fInit sh =                        -- negative of RHS of Poisson Equation = 0 for Laplace
    fromListUnboxed sh $ replicate (size sh) 1.0
                          

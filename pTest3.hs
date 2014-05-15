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
	  [stackSize,filename]	
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
     arrAns = exact shapeInit
  in
    do 
      (arrFinal, t) <- time $ fullMG h f bM bV
      putStr (prettyTime t)
      err :: Array U DIM2 Double 
          <- computeP 
             $ R.map abs
             $ R.zipWith (-) arrFinal arrAns
      putStrLn $ "max error = " Prelude.++ (show $ foldAllS max 0.0 err)
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
  = (Z :. n :. n) where n = 1+2^gridStackSize


{- Laplace Example with Exact solution
    
     http://www5.in.tum.de/lehre/vorlesungen/sci_comp/ws04/PDF/Lesson_09.pdf 
    
    (here x,y reversed and taking g_k = 0 except for particular k=1,2,... given by gCoeffs)
    
    >>>> Au = -f <<<< 

   u and f defined on the unit square (u is unknown)

   f(x,y) = 0 for Laplace
    
   Dirichlet Boundary Conditions

   u(x,0) = 0 and u(x,1) = 0            for 0<=x<=1
   u(0,y) = 0                           for 0<=y<=1
   u(1,y) = Sigma k=0 to infinity (g_k * sin(k*pi*y))           for 0<=y<=1


   Exact solution f
    
    u(x,y) = Sigma k=0 to infinity (sin (k*pi*y) * sinh(k*pi*x) * g_k / sinh(k*pi))

   Take g_k = 0 except for list of values for first few k (gCoeffs)
-}
   
distance:: Double -- length of sides for square area covered by grids
distance = 1.0 

gCoeffs :: [Double]  -- coefficients to play around with
gCoeffs = [5.0, 2.0, 3.0, 1.0]

gSigma :: Num b => (b -> Double) -> Double
gSigma f = sum [ f (fromIntegral i) * gCoeffs !! (i-1) | i <- [1..length gCoeffs]]


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
     [edgeRow,
      take ((m-2)*n) (cycle mainRow),
      edgeRow
     ]
     where (Z:.m:.n) = sh
           edgeRow = replicate n 0.0
           mainRow = 0.0: replicate (n-2) 1.0  Prelude.++ [0.0]

boundValues :: DIM2 -> Array U DIM2 Double        
boundValues sh =
  fromListUnboxed sh $ concat
     [
       take ((m-1)*n) (cycle mainRow),
       edgeRow
     ] 
     where (Z:.m:.n) = sh
           mainRow = replicate n 0.0
           edgeRow = Prelude.map (\j -> gSigma (\k -> sin(k*pi*j))) coords Prelude.++[0.0]
           coords = init (hCoords sh)

fInit :: DIM2 -> Array U DIM2 Double
fInit sh =                        -- negative of RHS of Poisson Equation = 0 for Laplace
    fromListUnboxed sh $ replicate (size sh) 0.0

exact :: DIM2 -> Array U DIM2 Double
exact sh =
      fromListUnboxed sh $
      concat $ Prelude.map row coords    
      where row i = Prelude.map (item i) coords
            item i j = gSigma (\k -> sin (k*pi*j) * sinh(k*pi*i) / sinh(k*pi))
            coords = hCoords sh
       
          
          
         


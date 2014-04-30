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



{- Example from p162 Iserles "A First Course in Numerical Analysis of Differential Equations"
   Values for poisson equation 
    
    >>>> Au = -f <<<< 

   u and f defined on the unit square (u is unknown)

   f(x,y) = -(x^2+y^2)

   Dirichlet Boundary Conditions

   u(x,0) = 0 and u(x,1) = 1/2 * x^2            for 0<=x<=1
   u(0,y) = sin (pi*x)                          for 0<=y<=1
   u(1,y) = e^(pi*x) * sin(pi*y) + 1/2 * y^2    for 0<=y<=1
-}

distance:: Double
distance = 1.0  --length of sides for square area covered by grids

fineGridSpaces:: Int -> Int
fineGridSpaces gridStackSize
  = 2^gridStackSize

fineGridShape:: Int -> DIM2
fineGridShape gridStackSize 
  = (Z :. n :. n) where n = 1+fineGridSpaces gridStackSize

writeHeatMapBMP :: String 
                -> Array U DIM2 Double
                -> IO()
writeHeatMapBMP filename arr
  = do let maxVal = foldAllS max 0.0 arr
           minVal = foldAllS min 0.0 arr
       arrImageOut <- computeP
                      $  R.map rgb8OfDouble
                      $  R.map (rampColorHotToCold minVal maxVal) arr
       writeImageToBMP filename arrImageOut

    
poissonTest :: Int -> String -> IO()
poissonTest gridStackSize filename =
 let intervals = fineGridSpaces gridStackSize

     shapeInit = fineGridShape gridStackSize

     hInit = distance / fromIntegral intervals
     
     boundMask = 
      fromListUnboxed shapeInit $ concat 
         [edgeRow,
          take ((intervals-1)*(intervals+1)) (cycle mainRow),
          edgeRow
         ]
         where edgeRow = replicate (intervals+1) 0.0
               mainRow = 0.0: replicate (intervals-1) 1.0  Prelude.++ [0.0]
     
     coordList = Prelude.map ((hInit*) . fromIntegral) [0..intervals] 
     
     boundValues =
      fromListUnboxed shapeInit $ concat
         [Prelude.map (\j -> sin (pi*j)) coordList, 
          concat $ Prelude.map mainRow $ tail $ init coordList,
          Prelude.map (\j -> exp pi * sin (pi*j) + 0.5*j^2) coordList
         ] 
         where mainRow i = replicate intervals 0.0 Prelude.++ [0.5*i^2]
     
     fInit =                        -- negative of RHS of Poisson Equation 
      fromListUnboxed shapeInit $ concat $ Prelude.map row coordList    
      where row i = Prelude.map (item i) coordList
            item i j = -(i^2 + j^2)

     exact =
      fromListUnboxed shapeInit $
      concat $ Prelude.map row coordList    
      where row i = Prelude.map (item i) coordList
            item i j = exp (pi*i) * sin (pi*j) + 0.5*i^2*j^2

     in
       do 
          (arrFinal, t) <- time $ fullMG hInit fInit boundMask boundValues
          putStr (prettyTime t)
          err :: Array U DIM2 Double 
              <- computeP 
                 $ R.map abs
                 $ R.zipWith (-) arrFinal exact
          putStrLn $ "max error = " Prelude.++ (show $ foldAllS max 0.0 err)
          writeHeatMapBMP filename arrFinal
          
          
          


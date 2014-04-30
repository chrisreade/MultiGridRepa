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


module MultiGrid_Poisson
       (multiGrid,fullMG)
 
where

import Data.Array.Repa				    as R
import Data.Array.Repa.Stencil			as R
import Data.Array.Repa.Stencil.Dim2		as R
import Data.Array.Repa.Eval             as R

import PoissonOps


steps1, steps2 :: Int
steps1 = 5       -- number of iterations for first step approximation operations
steps2 = 5       -- number of iterations for last step approximation operations

multiGrid :: Monad m =>
                  Double
                  -> Array U DIM2 Double
                  -> Array U DIM2 Double
                  -> Array U DIM2 Double
                  -> Array U DIM2 Double
                  -> m (Array U DIM2 Double)
                  
multiGrid !h !f !boundMask !boundValues !uInit 
 = if coarsest uInit
   then 
     do computeP $ approxOp h f boundMask boundValues uInit
   else
     do v  <- iterateSolver (approxOp h f boundMask boundValues) steps1 uInit
        r  <- computeP $ residualOp h f boundMask v
        boundMask' <- computeP $ coarsen boundMask
        r' <- computeP $ szipWith (*) boundMask' $ restrict r                            
        let zeros = zeroArray (extent r')
        err <- multiGrid (2*h) r' boundMask' zeros zeros
        vC  <- computeP $ szipWith (+) v  
                        $ szipWith (*) boundMask 
                        $ interpolate err
        iterateSolver (approxOp h f boundMask boundValues) steps2 vC 
             
fullMG :: Monad m =>
               Double
               -> Array U DIM2 Double
               -> Array U DIM2 Double
               -> Array U DIM2 Double
               -> m (Array U DIM2 Double)

fullMG !h !f !boundMask !boundValues
 = if coarsest boundValues
     then do computeP $ approxOp h f boundMask boundValues boundValues
     else do 
            f' <- computeP $ restrict f
            boundMask' <- computeP $ coarsen boundMask
            boundValues' <- computeP $ coarsen boundValues
            v'  <- fullMG (2*h) f' boundMask' boundValues'  
            v <- computeP $ szipWith (+) boundValues  
                          $ szipWith (*) boundMask 
                          $ interpolate v'
            multiGrid h f boundMask boundValues v


{-# INLINE coarsest #-}
coarsest :: Array U DIM2 Double -> Bool
coarsest !a = m<4 where (Z :. _ :. m) = extent a

{-# INLINE zeroArray #-}
zeroArray  :: Shape sh => sh -> Array U sh Double
zeroArray !sh = fromListUnboxed sh $ replicate (size sh) 0.0

{-# INLINE iterateSolver #-}
iterateSolver :: (Monad m, Source r2 e, Load r1 sh e, Target r2 e) =>
                 (Array r2 sh e -> Array r1 sh e)
                 -> Int 
                 -> Array r2 sh e 
                 -> m (Array r2 sh e)
iterateSolver !opA !steps !arrInit
 = go steps arrInit
  where
    go 0 !arr = return arr
    go n !arr 
         = do   arr' <- computeP $ opA arr
                go (n - 1) arr'

restStencil :: Stencil DIM2 Double
restStencil =   [stencil2|  1 2 1         
                            2 4 2 
                            1 2 1 |]

sum4Stencil :: Stencil DIM2 Double    
sum4Stencil = [stencil2|  0 0 0         
                          0 1 1 
                          0 1 1 |] 

{-# INLINE restrict #-}       
restrict :: Array U DIM2 Double -> Array D DIM2 Double
restrict !arr 
  | odd n && odd m
    = coarsen
      $ smap (/16)  
      $ mapStencil2 BoundClamp restStencil arr
  | otherwise
    = coarsen
      $ smap (/4)  
      $ mapStencil2 BoundClamp sum4Stencil arr
  where _ :. m :. n = extent arr

{-# INLINE coarsen #-}        
coarsen :: Source r a => Array r DIM2 a -> Array D DIM2 a
coarsen !arr = traverse arr   -- i+1 and j+1 to deal with odd extents correctly
                        (\ (e :. i :. j) -> (e :. (i+1) `div` 2 :. (j+1) `div` 2))
                        (\get (e :. i :. j) -> get (e :. 2*i :. 2*j))
                         

{-# INLINE inject4 #-} 
inject4 :: Source r a => Array r DIM2 a -> Array D DIM2 a
inject4 !arr = traverse arr -- mod 2s deal with odd extents
                        (\ (e :. i :. j) -> (e :. 2*i - (i `mod` 2) :. 2*j - (j `mod` 2)))
                        (\get (e :. i :. j) -> get(e :. i `div` 2 :. j `div` 2))

{-# INLINE interpolate #-} 
interpolate :: Array U DIM2 Double -> Array D DIM2 Double
interpolate !arr 
   | odd m && odd n
       = traverse arr 
                  (\ (e :. i :. j) -> (e :. 2*i - (i `mod` 2) :. 2*j - (j `mod` 2)))
                  (\get (e :. i :. j) -> case () of
                   _ | even i && even j     -> get(e :. i `div` 2 :. j `div` 2)
                     | even i               -> (0.5)*(get(e :. i `div` 2 :. j `div` 2)
                                                      + get(e :. i `div` 2 :. (j `div` 2)+1)) 
                      -- odd i
                     | even j               -> (0.5)*(get(e :. i `div` 2 :. j `div` 2)
                                                      + get(e :. (i `div` 2)+1 :. j `div` 2)) 
                      -- odd i and j
                     | otherwise            -> (0.25)*(get(e :. i `div` 2 :. j `div` 2)
                                                      + get(e :. i `div` 2 :. (j `div` 2)+1)
                                                      + get(e :. (i `div` 2)+1 :. j `div` 2)
                                                      + get(e :. (i `div` 2)+1 :. (j `div` 2)+1))
                    )
   | otherwise = inject4 arr
  where _ :. n :. m = extent arr

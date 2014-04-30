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

module PoissonOps

where

import Data.Array.Repa				    as R
import Data.Array.Repa.Stencil			as R
import Data.Array.Repa.Stencil.Dim2		as R

{-# INLINE approxOp #-} 
approxOp :: Double
            -> Array U DIM2 Double
            -> Array U DIM2 Double
            -> Array U DIM2 Double
            -> Array U DIM2 Double
            -> Array (TR PC5) DIM2 Double               

approxOp !h !f !boundMask !boundValues !arr
       =   szipWith (+) boundValues
         $ szipWith (*) boundMask
         $ smap (/4)
         $ szipWith (+) (R.map (*hSq) f )
         $ mapStencil2 (BoundConst 0) 
            [stencil2|   0 1 0
                         1 0 1 
                         0 1 0 |] arr
         where hSq = h*h

{-# INLINE residualOp #-}
residualOp :: Double
              -> Array U DIM2 Double
              -> Array U DIM2 Double
              -> Array U DIM2 Double
              -> Array (TR PC5) DIM2 Double

residualOp !h !f !boundMask !v
       =   szipWith (*) boundMask
         $ szipWith (+) f
         $ smap (*hFactor)
         $ mapStencil2 (BoundConst 0) 
            [stencil2|   0 1  0
                         1 -4 1 
                         0 1  0 |] v
         where hFactor = 1/(h*h)

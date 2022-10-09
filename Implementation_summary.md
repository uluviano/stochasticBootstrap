Implementation summary:
   ```
   Parameter       & Value(s)   & Description\\
  \texttt{ nZ                }& $200     $  & Number of points taken in the $z-$plane.\\
  \texttt{ lambda0           }& $0.42    $  & Size of the region in which the points are sampled.\\
  \texttt{ Dstar             }& $7-15$     & $\Delta^*$ used in the normalization of the action \eqref{}\\
  \texttt{ NT                }& $ 2 \times 10^8$ &  Number of steps (MC time). \\
  \texttt{ positiveOPEs      }& \texttt{true }&     Enforcement of real  OPE coefficients. \\
  \texttt{ externalInOPE     }& \texttt{false} &    Fixing the first exchanged scalar to be equal to the external operator.\\
  \texttt{ Temp              }& $ 0.3-0.9  $ &    Temperature ($T$) of the Metropolis Monte Carlo algorithm.\\
  \texttt{ MCstep            }& $ 0.001    $ &    Overall scale factor for the MC-step.   \\
  \texttt{ wall              }& $ 10000.   $ &    Constant for the quadratic penalty that enforces the bounds on $\Delta_i$. This guarantees analiticity of the potential. \\
  \texttt{ invSampleDensity  } & $ 1000     $ &  The MC saves one out of   \texttt{ invSampleDensity  } frames.\\
  \texttt{ seedInput         }& $ 211      $ & The seed of the pRNG. This guarantees reproducibility.\\
  \texttt{ nops              }& $ 3 - 20   $ &  Number of operators (besides the identity).\\
   ```
			
			 
			 
			     
			     
			  
			      
			    
			   
/*
 * DTW.java   Jul 14, 2004
 *
 * Copyright (c) 2004 Stan Salvador
 * stansalvador@hotmail.com
 */

/*
 * This class contains a "light" implementation of DTW. There are no motion weights, and you the main option is whether you want DTW to be flexible or not.
 *
 */
package com.dtw;

import com.timeseries.TimeSeries;
import com.util.DistanceFunction;
import com.matrix.ColMajorCell;

import java.util.Arrays;
import java.util.Iterator;


public class LightDTW
{

   // FUNCTIONS
   public static double calcWarpCost(WarpPath path, TimeSeries tsI, TimeSeries tsJ, DistanceFunction distFn)
   {
      double totalCost = 0.0;

      for (int p=0; p<path.size(); p++)
      {
         final com.matrix.ColMajorCell currWarp = path.get(p);
         totalCost += distFn.calcDistance(tsI.getMeasurementVector(currWarp.getCol()),
                                          tsJ.getMeasurementVector(currWarp.getRow()));
      }

      return totalCost;
   }


   // Dynamic Time Warping where the warp path is not needed, an alternate implementation can be used that does not
   //    require the entire cost matrix to be filled and only needs 2 columns to be stored at any one time.
 /*  public static double getWarpDistBetween(TimeSeries tsI, TimeSeries tsJ, DistanceFunction distFn)
   {
      // The space complexity is 2*tsJ.size().  Dynamic time warping is symmetric so switching the two time series
      //    parameters does not effect the final warp cost but can reduce the space complexity by allowing tsJ to
      //    be set as the shorter time series and only requiring 2 columns of size |tsJ| rather than 2 larger columns of
      //    size |tsI|.
      if (tsI.size() < tsJ.size())
         return getWarpDistBetween(tsJ, tsI, distFn);


      double[] lastCol = new double[tsJ.size()];
      double[] currCol = new double[tsJ.size()];
      final int maxI = tsI.size()-1;
      final int maxJ = tsJ.size()-1;

      // Calculate the values for the first column, from the bottom up.
      currCol[0] = distFn.calcDistance(tsI.getMeasurementVector(0), tsJ.getMeasurementVector(0));  // first cell
      for (int j=1; j<=maxJ; j++)  // the rest of the first column
         currCol[j] = currCol[j-1] + distFn.calcDistance(tsI.getMeasurementVector(0), tsJ.getMeasurementVector(j));

      for (int i=1; i<=maxI; i++)   // i = columns
      {
         // Swap the references between the two arrays.
         final double[] temp = lastCol;
         lastCol = currCol;
         currCol = temp;

         // Calculate the value for the bottom row of the current column
         //    (i,0) = LocalCost(i,0) + GlobalCost(i-1,0)
         currCol[0] = lastCol[0] + distFn.calcDistance(tsI.getMeasurementVector(i), tsJ.getMeasurementVector(0));

         for (int j=1; j<=maxJ; j++)  // j = rows
         {
            // (i,j) = LocalCost(i,j) + minGlobalCost{(i-1,j),(i-1,j-1),(i,j-1)}
            final double minGlobalCost = Math.min(lastCol[j], Math.min(lastCol[j-1], currCol[j-1]));
            currCol[j] = minGlobalCost + distFn.calcDistance(tsI.getMeasurementVector(i), tsJ.getMeasurementVector(j));
         }  // end for loop
      }  // end for loop

      // Minimum Cost is at (maxI,maxJ)
      return currCol[maxJ];
   }  // end getWarpDistBetween(..) */


   public static WarpPath getWarpPathBetween(TimeSeries tsI, TimeSeries tsJ, DistanceFunction distFn, boolean flexible)
   {
      return DynamicTimeWarp(tsI, tsJ, distFn, flexible).getPath();
   }


   public static TimeWarpInfo getWarpInfoBetween(TimeSeries tsI, TimeSeries tsJ, DistanceFunction distFn, boolean flexible)
   {
      return DynamicTimeWarp(tsI, tsJ, distFn, flexible);
   }


   private static TimeWarpInfo DynamicTimeWarp(TimeSeries tsI, TimeSeries tsJ, DistanceFunction distFn, boolean flexible)
   {
      //     COST MATRIX:
      //   5|_|_|_|_|_|_|E| E = min Global Cost
      //   4|_|_|_|_|_|_|_| S = Start point
      //   3|_|_|_|_|_|_|_| each cell = min global cost to get to that point
      // j 2|_|_|_|_|_|_|_|
      //   1|_|_|_|_|_|_|_|
      //   0|S|_|_|_|_|_|_|
      //     0 1 2 3 4 5 6
      //            i
      //   access is M(i,j)... column-row

      final double[][] costMatrix = new double[tsI.size()][tsJ.size()];
      final int maxI = tsI.size()-1;
      final int maxJ = tsJ.size()-1;

      // Fill the values in the for i = 0. If DTW is flexible, just put the distance of tsI[0] and tsJ[j], because 
      // flexible DTW ends warping whenever i = 0.

      costMatrix[0][0] = distFn.calcDistance(tsI.getMeasurementVector(0),
                                       tsJ.getMeasurementVector(0));

      if (flexible) {
	      for (int j=1; j<=maxJ; j++)
		 costMatrix[0][j] = distFn.calcDistance(tsI.getMeasurementVector(0),
		                                                             tsJ.getMeasurementVector(j));
      } else {
	      for (int j=1; j<=maxJ; j++)
		 costMatrix[0][j] = costMatrix[0][j-1] + distFn.calcDistance(tsI.getMeasurementVector(0),
		                                                             tsJ.getMeasurementVector(j));
      }

      for (int i=1; i<=maxI; i++)
      {
         // Calculate the value for the bottom row of the current column
         //    (i,0) = LocalCost(i,0) + GlobalCost(i-1,0).
         // If j = 0 and i > 0, you can only go to i-1, j.
         costMatrix[i][0] = costMatrix[i-1][0] + distFn.calcDistance(tsI.getMeasurementVector(i),
                                                                     tsJ.getMeasurementVector(0));

         for (int j=1; j<=maxJ; j++)  // j = rows
         {
            // (i,j) = LocalCost(i,j) + minGlobalCost{(i-1,j),(i-1,j-1),(i,j-1)}
            final double minGlobalCost = Math.min(costMatrix[i-1][j],
                                                  Math.min(costMatrix[i-1][j-1],
                                                           costMatrix[i][j-1]));
            costMatrix[i][j] = minGlobalCost + distFn.calcDistance(tsI.getMeasurementVector(i),
                                                                   tsJ.getMeasurementVector(j));
         }  // end for loop
      }  // end for loop

      // Minimum Cost is at (maxI,maxJ) for Inflexible, and (maxI, bestJ), 0 <= bestJ <= maxJ, for Flexible.
      double minimumCost = Double.POSITIVE_INFINITY;
      int bestJ = -1;
      if (flexible) {
	for (int j = 0; j <= maxJ; j++) 
		if (costMatrix[maxI][j] < minimumCost)  {
			minimumCost = costMatrix[maxI][j];
			bestJ = j;
		}

      } else {
	 minimumCost = costMatrix[maxI][maxJ];
         bestJ = maxJ;
      }


      // Find the Warp Path by searching the matrix from the solution at
      //    (maxI, bestJ) to the beginning at (0,starting_j).  At each step move through
      //    the matrix 1 step left, down, or diagonal, whichever has the
      //    smallest cost.  Favour diagonal moves and moves towards the i==j
      //    axis to break ties.
      final WarpPath minCostPath = new WarpPath(maxI+maxJ-1);
      int i = maxI;
      int j = bestJ;
      minCostPath.addFirst(i, j);
     
      // If DTW is Flexible, end your backtracking as soon as i == 0; (0, j) will be the start of your warp path.
      if (flexible) {
	      while (i>0)
	      {
	         // Find the costs of moving in all three possible directions (left,
		 //    down, and diagonal (down and left at the same time).
		 final double diagCost;
		 final double leftCost;
		 final double downCost;

                 leftCost = costMatrix[i-1][j];

		 if (j > 0) {
                    downCost = costMatrix[i][j-1];
		    diagCost = costMatrix[i-1][j-1];
		 } else {
                    downCost = Double.POSITIVE_INFINITY;
		    diagCost = Double.POSITIVE_INFINITY;
                 }

		 // Determine which direction to move in.  Prefer moving diagonally and
		 //    moving towards the i==j axis of the matrix if there are ties.
		 if ((diagCost<=leftCost) && (diagCost<=downCost))
		 {
		    i--;
		    j--;
		 }
		 else if ((leftCost<diagCost) && (leftCost<downCost))
		    i--;
		 else if ((downCost<diagCost) && (downCost<leftCost))
		    j--;
		 else if (i <= j)  // leftCost==rightCost > diagCost
		    j--;
		 else   // leftCost==rightCost > diagCost
		    i--;

		 // Add the current step to the warp path.
		 minCostPath.addFirst(i, j);
	      }  // end while loop
      } else {
	      while ((i>0) || (j>0))
	      {
		 // Find the costs of moving in all three possible directions (left,
		 //    down, and diagonal (down and left at the same time).
		 final double diagCost;
		 final double leftCost;
		 final double downCost;

		 if ((i>0) && (j>0))
		    diagCost = costMatrix[i-1][j-1];
		 else
		    diagCost = Double.POSITIVE_INFINITY;

		 if (i > 0)
		    leftCost = costMatrix[i-1][j];
		 else
		    leftCost = Double.POSITIVE_INFINITY;

		 if (j > 0)
		    downCost = costMatrix[i][j-1];
		 else
		    downCost = Double.POSITIVE_INFINITY;

		 // Determine which direction to move in.  Prefer moving diagonally and
		 //    moving towards the i==j axis of the matrix if there are ties.
		 if ((diagCost<=leftCost) && (diagCost<=downCost))
		 {
		    i--;
		    j--;
		 }
		 else if ((leftCost<diagCost) && (leftCost<downCost))
		    i--;
		 else if ((downCost<diagCost) && (downCost<leftCost))
		    j--;
		 else if (i <= j)  // leftCost==rightCost > diagCost
		    j--;
		 else   // leftCost==rightCost > diagCost
		    i--;

		 // Add the current step to the warp path.
		 minCostPath.addFirst(i, j);
	     }  // end while loop
      }

      return new TimeWarpInfo(minimumCost, minCostPath);
   }  // end DynamicTimeWarp(..)


   /*
   public static double getWarpDistBetween(TimeSeries tsI, TimeSeries tsJ, SearchWindow window, DistanceFunction distFn)
   {
      //     COST MATRIX:
      //   5|_|_|_|_|_|_|E| E = min Global Cost
      //   4|_|_|_|_|_|_|_| S = Start point
      //   3|_|_|_|_|_|_|_| each cell = min global cost to get to that point
      // j 2|_|_|_|_|_|_|_|
      //   1|_|_|_|_|_|_|_|
      //   0|S|_|_|_|_|_|_|
      //     0 1 2 3 4 5 6
      //            i
      //   access is M(i,j)... column-row
      final CostMatrix costMatrix = new PartialWindowMatrix(window);
      final int maxI = tsI.size()-1;
      final int maxJ = tsJ.size()-1;

      // Get an iterator that traverses the window cells in the order that the cost matrix is filled.
      //    (first to last row (1..maxI), bottom to top (1..MaxJ)
      final Iterator matrixIterator = window.iterator();

      while (matrixIterator.hasNext())
      {
         final ColMajorCell currentCell = (ColMajorCell)matrixIterator.next();  // current cell being filled
         final int i = currentCell.getCol();
         final int j = currentCell.getRow();

         if ( (i==0) && (j==0) )      // bottom left cell (first row AND first column)
            costMatrix.put(i, j, distFn.calcDistance(tsI.getMeasurementVector(0), tsJ.getMeasurementVector(0)));
         else if (i == 0)             // first column
         {
            costMatrix.put(i, j, distFn.calcDistance(tsI.getMeasurementVector(0), tsJ.getMeasurementVector(j)) +
                                 costMatrix.get(i, j-1));
         }
         else if (j == 0)             // first row
         {
            costMatrix.put(i, j, distFn.calcDistance(tsI.getMeasurementVector(i), tsJ.getMeasurementVector(0)) +
                                 costMatrix.get(i-1, j));
         }
         else                         // not first column or first row
         {
            final double minGlobalCost = Math.min(costMatrix.get(i-1, j),
                                                  Math.min(costMatrix.get(i-1, j-1),
                                                           costMatrix.get(i, j-1)));
            costMatrix.put(i, j, minGlobalCost + distFn.calcDistance(tsI.getMeasurementVector(i),
                                                                     tsJ.getMeasurementVector(j)));
         }  // end if
      }  // end while loop

      // Minimum Cost is at (maxI, maxJ)
      return costMatrix.get(maxI, maxJ);

   }  // end getWarpDistBetween(...)
   */

   public static WarpPath getWarpPathBetween(TimeSeries tsI, TimeSeries tsJ, SearchWindow window, DistanceFunction distFn, boolean flexible)
   {
      return constrainedTimeWarp(tsI, tsJ, window, distFn, flexible).getPath();
   }


   public static TimeWarpInfo getWarpInfoBetween(TimeSeries tsI, TimeSeries tsJ, SearchWindow window, DistanceFunction distFn, boolean flexible)
   {
      return constrainedTimeWarp(tsI, tsJ, window, distFn, flexible);
   }


   private static TimeWarpInfo constrainedTimeWarp(TimeSeries tsI, TimeSeries tsJ, SearchWindow window, DistanceFunction distFn, boolean flexible)
   {
      //     COST MATRIX:
      //   5|_|_|_|_|_|_|E| E = min Global Cost
      //   4|_|_|_|_|_|_|_| S = Start point
      //   3|_|_|_|_|_|_|_| each cell = min global cost to get to that point
      // j 2|_|_|_|_|_|_|_|
      //   1|_|_|_|_|_|_|_|
      //   0|S|_|_|_|_|_|_|
      //     0 1 2 3 4 5 6
      //            i
      //   access is M(i,j)... column-row

      final WindowMatrix costMatrix = new WindowMatrix(window);
      final int maxI = tsI.size()-1;
      final int maxJ = tsJ.size()-1;

      // Get an iterator that traverses the window cells in the order that the cost matrix is filled.
      //    (first to last row (1..maxI), bottom to top (1..MaxJ)
      final Iterator matrixIterator = window.iterator();

      if (flexible) { 
	    while (matrixIterator.hasNext())
	      {
		 final ColMajorCell currentCell = (ColMajorCell)matrixIterator.next();  // current cell being filled
		 final int i = currentCell.getCol();
		 final int j = currentCell.getRow();

		
		 if (i == 0)      // If i = 0, you don't have to move anywhere.
		 {
		    costMatrix.put(i, j, distFn.calcDistance(tsI.getMeasurementVector(0), tsJ.getMeasurementVector(j)));
		 }
		 else if (j == 0) // If j = 0, you can only go in the direction i-1, j.
		 {
		    costMatrix.put(i, j, distFn.calcDistance(tsI.getMeasurementVector(i), tsJ.getMeasurementVector(0)) +
		                         costMatrix.get(i-1, j));
		 }
		 else             // If i > 0 and j > 0, you can go in any of the three directions. 
		 {
		    final double minGlobalCost = Math.min(costMatrix.get(i-1, j),
		                                          Math.min(costMatrix.get(i-1, j-1),
		                                                   costMatrix.get(i, j-1)));
		    costMatrix.put(i, j, minGlobalCost + distFn.calcDistance(tsI.getMeasurementVector(i),
		                                                             tsJ.getMeasurementVector(j)));
		 }  // end if
	      }  // end while loop
      } else {
	      while (matrixIterator.hasNext())
	      {
		 final ColMajorCell currentCell = (ColMajorCell)matrixIterator.next();  // current cell being filled
		 final int i = currentCell.getCol();
		 final int j = currentCell.getRow();

		 if ( (i==0) && (j==0) )      // The start of the warp path
		    costMatrix.put(i, j, distFn.calcDistance(tsI.getMeasurementVector(0), tsJ.getMeasurementVector(0)));

		 else if (i == 0)             // Can only move to i, j-1
		 {
		    costMatrix.put(i, j, distFn.calcDistance(tsI.getMeasurementVector(0), tsJ.getMeasurementVector(j)) +
		                         costMatrix.get(i, j-1));
		 }
		 else if (j == 0)             // Can only move to i-1, j
		 {
		    costMatrix.put(i, j, distFn.calcDistance(tsI.getMeasurementVector(i), tsJ.getMeasurementVector(0)) +
		                         costMatrix.get(i-1, j));
		 }
		 else                         // Can move in any of the 3 directions
		 {
		    final double minGlobalCost = Math.min(costMatrix.get(i-1, j),
		                                          Math.min(costMatrix.get(i-1, j-1),
		                                                   costMatrix.get(i, j-1)));
		    costMatrix.put(i, j, minGlobalCost + distFn.calcDistance(tsI.getMeasurementVector(i),
		                                                             tsJ.getMeasurementVector(j)));
		 }  // end if
	      }  // end while loop
      }

      // Minimum Cost is at (maxI, bestJ), bestJ = maxJ for Inflexible, 0 <= bestJ <= maxJforI(maxI) for Flexible
      // maxJforI(maxI) is the greatest value of J on the last column of the window.
      double minimumCost = Double.POSITIVE_INFINITY;
      int bestJ = -1;
      if (flexible) {
	  for (int j = window.minJforI(maxI); j <= window.maxJforI(maxI); j++)
		if (costMatrix.get(maxI, j) < minimumCost) {
			minimumCost = costMatrix.get(maxI, j);
			bestJ = j;
		}
      } else {     
	  minimumCost = costMatrix.get(maxI, maxJ);
          bestJ = maxJ; // For Inflexible DTW, tsJ is completely used, unlike for Flexible DTW. 
      }


      // Find the Warp Path by searching the matrix from the solution at
      //    (maxI, bestJ) to the beginning at (0,starting_j).  At each step move through
      //    the matrix 1 step left, down, or diagonal, whichever has the
      //    smallest cost.  Favoer diagonal moves and moves towards the i==j
      //    axis to break ties.
      final WarpPath minCostPath = new WarpPath(maxI+maxJ-1);
      int i = maxI;
      int j = bestJ;
      minCostPath.addFirst(i, j);

      if (flexible) {
	while (i>0)
      	{
		 // Find the costs of moving in all three possible directions (left,
		 //    down, and diagonal (down and left at the same time).
		 final double diagCost;
		 final double leftCost;
		 final double downCost;

		 leftCost = costMatrix.get(i-1, j);

		 if (j>0) {
                    downCost = costMatrix.get(i, j-1);
		    diagCost = costMatrix.get(i-1, j-1);
		 } else {
                    downCost = Double.POSITIVE_INFINITY;
		    diagCost = Double.POSITIVE_INFINITY;
	         }

		 // Determine which direction to move in.  Prefer moving diagonally and
		 //    moving towards the i==j axis of the matrix if there are ties.
		 if ((diagCost<=leftCost) && (diagCost<=downCost))
		 {
		    i--;
		    j--;
		 }
		 else if ((leftCost<diagCost) && (leftCost<downCost))
		    i--;
		 else if ((downCost<diagCost) && (downCost<leftCost))
		    j--;
		 else if (i <= j)  // leftCost==rightCost > diagCost
		    j--;
		 else   // leftCost==rightCost > diagCost
		    i--;

		 // Add the current step to the warp path.
		 minCostPath.addFirst(i, j);
    	}  // end while loop

      } else {
	      while ((i>0) || (j>0))
	      {
		 // Find the costs of moving in all three possible directions (left,
		 //    down, and diagonal (down and left at the same time).
		 final double diagCost;
		 final double leftCost;
		 final double downCost;

		 if ((i>0) && (j>0))
		    diagCost = costMatrix.get(i-1, j-1);
		 else
		    diagCost = Double.POSITIVE_INFINITY;

		 if (i > 0)
		    leftCost = costMatrix.get(i-1, j);
		 else
		    leftCost = Double.POSITIVE_INFINITY;

		 if (j > 0)
		    downCost = costMatrix.get(i, j-1);
		 else
		    downCost = Double.POSITIVE_INFINITY;

		 // Determine which direction to move in.  Prefer moving diagonally and
		 //    moving towards the i==j axis of the matrix if there are ties.
		 if ((diagCost<=leftCost) && (diagCost<=downCost))
		 {
		    i--;
		    j--;
		 }
		 else if ((leftCost<diagCost) && (leftCost<downCost))
		    i--;
		 else if ((downCost<diagCost) && (downCost<leftCost))
		    j--;
		 else if (i <= j)  // leftCost==rightCost > diagCost
		    j--;
		 else   // leftCost==rightCost > diagCost
		    i--;

		 // Add the current step to the warp path.
		 minCostPath.addFirst(i, j);
	      }  // end while loop
      }

      // Free any rescources associated with the costMatrix (a swap file may have been created if the swa file did not
      //    fit into main memory).
      costMatrix.freeMem();

      return new TimeWarpInfo(minimumCost, minCostPath);
   }  // end ConstrainedTimeWarp

}  // end class DtwTest

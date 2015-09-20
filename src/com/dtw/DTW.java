/*
 * DTW.java   Jul 14, 2004
 *
 * Copyright (c) 2004 Stan Salvador
 * stansalvador@hotmail.com
 */

package com.dtw;

import com.timeseries.TimeSeries;
import com.util.DistanceFunction;
import com.util.EuclideanDistance;
import com.util.DotProductDistance;
import com.matrix.ColMajorCell;

import com.util.DistanceFunctionFactory;
import java.util.Arrays;
import java.util.Iterator;
import java.lang.Math;

import java.lang.String;

public class DTW
{
   /*
   //This does not return the right cost for weighted DTW b/c it doesn't include the weights
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
      return -1;
   }*/


  /* PARAMETERS 

	 - tsI, tsJ are the input TimeSeries to warp.

	 - distFn is the DistanceFunction you want to use. See src/com/util for the DistanceFcuntion classes.

	 - vcost is a constant cost associated with moving from (i, j) to (i-1, j), duplicating a database frame.

	 - hcost is a constant cost associated with moving from (i, j) to (i, j-1), duplicating a query frame. 

	 - motionWeight is a multiplicative factor for the motion weights, which are penalties calculated during warping, based on query and database dynamics.

	 - motionDist is either "Dot" or "Euc", indicating the distance measure to use for calculating the motion weights.

	 - normalize indicates whether to normalize the warp path score by the warp path length. This is useful in flexible matching, where all
	   the query frames must be used, but not all the database frames, to avoid the inherent bias against horizontal movements, which would increase the warp path length.
           The algorithm doesn't compute the optimal normalized path.

  */

   public static WarpPath getWarpPathBetween(TimeSeries tsI, TimeSeries tsJ, DistanceFunction distFn)
   {
      return DynamicTimeWarp(tsI, tsJ, distFn, 0, 0, 0, "Dot", false).getPath();
   }

   public static WarpPath getWarpPathBetween(TimeSeries tsI, TimeSeries tsJ, DistanceFunction distFn, double vcost, double hcost, double motionWeight, String motionDist, boolean normalize)
   {
      return DynamicTimeWarp(tsI, tsJ, distFn, vcost, hcost, motionWeight, motionDist, normalize).getPath();
   }

   public static TimeWarpInfo getWarpInfoBetween(TimeSeries tsI, TimeSeries tsJ, DistanceFunction distFn)
   {
      return DynamicTimeWarp(tsI, tsJ, distFn, 0, 0, 0, "Dot", false);
   }

   public static TimeWarpInfo getWarpInfoBetween(TimeSeries tsI, TimeSeries tsJ, DistanceFunction distFn, double vcost, double hcost, double motionWeight, String motionDist, boolean normalize)
   {
      return DynamicTimeWarp(tsI, tsJ, distFn, vcost, hcost, motionWeight, motionDist, normalize);
   }


   private static TimeWarpInfo DynamicTimeWarp(TimeSeries tsI, TimeSeries tsJ, DistanceFunction distFn, double vcost, double hcost, double motionWeight, String motionDist, boolean normalize)
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

//double s;
//double e;

//s = System.currentTimeMillis();
      double[][] costMatrix = null;
      double[][] normCostMatrix = null;
      int[][] pathLength = null;

      // For normalization, we need to keep track of normalized cost and path length, instead of total cost.
      if (normalize) {
	 normCostMatrix = new double[tsI.size()][tsJ.size()];
	 pathLength = new int[tsI.size()][tsJ.size()];
      }
      else
	 costMatrix = new double[tsI.size()][tsJ.size()];

//e = System.currentTimeMillis();
//System.out.println("Allocation: " + (e-s));


      final int maxI = tsI.size()-1;
      final int maxJ = tsJ.size()-1;
      int dim = tsI.getMeasurementVector(0).length;
      // Create a zero vector, for calculating Euclidean distance motion weights when duplicating a frame (i,j) -> (i-1, j) or (i, j-1).
      double[] zero = null;
      if (motionDist.equals("Euc")) {
	zero = new double[dim];
        for (int z = 0; z < dim; z++)
          zero[z] = 0;
      }

      // These are used for calculating distances for the motion weights.
      DistanceFunction mot = DistanceFunctionFactory.getDistFnByName(motionDist);

//s = System.currentTimeMillis();
      // Precompute the difference between each frame and the frame before; these values are used in the motion weights.  
      double[][] tsIDiffs = new double[maxI][dim];
      double[][] tsJDiffs = new double[maxJ + 1][dim];	
      for (int tsind = 0; tsind < maxI; tsind++)
          tsIDiffs[tsind] = arrayMinus(tsI.getMeasurementVector(tsind+1), tsI.getMeasurementVector(tsind));
      for (int tsjind = 0; tsjind <= maxJ; tsjind++) {
          if (tsjind != 0)
	      tsJDiffs[tsjind] = arrayMinus(tsJ.getMeasurementVector(tsjind), tsJ.getMeasurementVector(tsjind-1));
	  else
	      tsJDiffs[tsjind] = tsJ.getMeasurementVector(tsjind); // If j == 0, consider the frame before J's 0th frame to be the zero vector.
      }
//e = System.currentTimeMillis();
//System.out.println("Velocity: " + (e-s));


//s = System.currentTimeMillis();
      // For cells (0, j), fill with the "difference score" between the first query frame and database frame j.
      // A warp path must end at one of these cells, because all query frames must be used, but not all database frames. 
      if (normalize) {
	 for (int j = 0; j <= maxJ; j++) {
		normCostMatrix[0][j] = distFn.calcDistance(tsI.getMeasurementVector(0),
                                                                     tsJ.getMeasurementVector(j));
		pathLength[0][j] = 1;
	}
      }	else {
	 for (int j = 0; j <= maxJ; j++)
		costMatrix[0][j] = distFn.calcDistance(tsI.getMeasurementVector(0),
                                                                     tsJ.getMeasurementVector(j));
      }
//e = System.currentTimeMillis();
//System.out.println("First row: " + (e-s));

//s = System.currentTimeMillis();
      // Here, the cost matrix is filled, using bottom-up Dynamic Programming.
      for (int i=1; i<=maxI; i++)
      {
         // First, consider j = 0. The only transition that can be made is (i, j) -> (i-1, j). Thus, cells (i, 0) on the matrix can be filled using only the value from cell (i-1, 0).
	 
	 // Calculate iCost, the *motion weight* cost for going from (i, j) -> (i-1, j).
         double[] difI = tsIDiffs[i-1];
	 double iCost = 0;
	 if (motionDist.equals("Euc"))
	 	iCost = motionWeight*mot.calcDistance(difI, zero);
	 else if (motionDist.equals("Dot")) {
		double[] difJ;
		difJ = tsJ.getMeasurementVector(0); //for j = 0, consider the frame before J's 0th frame to be the zero vector
		double[] jEpsilon = normalizeVector(momentumEpsilon(difJ));		
		iCost = motionWeight*mot.calcDistance(jEpsilon, difI);
	 } else
	 	iCost = 0;

	 // Fill in the cost matrix cell (i, 0).
	 if (normalize) {
		normCostMatrix[i][0] = (normCostMatrix[i-1][0]*pathLength[i-1][0] + motionWeight*distFn.calcDistance(tsI.getMeasurementVector(i), tsJ.getMeasurementVector(0)) + vcost + iCost ) / //!!!
					(pathLength[i-1][0] + 1);
		pathLength[i][0] = pathLength[i-1][0] + 1;
	 } else {
        	 costMatrix[i][0] = costMatrix[i-1][0] + distFn.calcDistance(tsI.getMeasurementVector(i),
                                                                     tsJ.getMeasurementVector(0)) + vcost + iCost;
	 }

	 // Consider j > 0.
         for (int j=1; j<=maxJ; j++)
         {
	    final double[] difJ = tsJDiffs[j-1];
	    
	    final double diagCost;
	    final double jCost;

	    // Re-calculate iCost for different values of jEpsilon.
	    if (motionDist.equals("Dot")) {
		double[] jEpsilon = normalizeVector(momentumEpsilon(difJ));		
		iCost = motionWeight*mot.calcDistance(jEpsilon, difI);
	    }

	    // Calculate the *motion weight* cost for transitions (i, j) -> (i-1, j-1) [diagCost] and (i, j) -> (i, j-1) [jCost] .
	    if (motionDist.equals("Euc"))
	        diagCost = motionWeight*mot.calcDistance(difI, difJ);
	    else if (motionDist.equals("Dot"))
		diagCost = 1*motionWeight*mot.calcDistance(normalizeVector(difI), normalizeVector(difJ));
	    else
		diagCost = 0;

	    if (motionDist.equals("Euc"))
		jCost = motionWeight*mot.calcDistance(difJ, zero);
	    else if (motionDist.equals("Dot")) {	
		double[] iEpsilon = normalizeVector(momentumEpsilon(difI));		
		jCost = motionWeight*mot.calcDistance(difJ, iEpsilon);	
	    } else
		jCost = 0;
	
  	    // Calculate the *full costs* for each of the three possible transitions, and choose the one that minimizes the cost of the path from cell (i,j) to a cell (0, j).  
	    if (normalize) {

		double myDist = distFn.calcDistance(tsI.getMeasurementVector(i), tsJ.getMeasurementVector(j));
		double fullDiagCost = (pathLength[i-1][j-1]*normCostMatrix[i-1][j-1] + diagCost + myDist) / (pathLength[i-1][j-1] + 1); 
		double fullICost = (pathLength[i-1][j]*normCostMatrix[i-1][j] + vcost + iCost + motionWeight*myDist) / (pathLength[i-1][j] + 1); //!!!
		double fullJCost = (pathLength[i][j-1]*normCostMatrix[i][j-1] + hcost + jCost + motionWeight*myDist) / (pathLength[i][j-1] + 1); //!!!
		
		final double minGlobalCost = Math.min(fullDiagCost, Math.min(fullICost, fullJCost));
		if (Double.compare(minGlobalCost + 0.0, fullDiagCost + 0.0) == 0) // +0.0 is to avoid a case where -0.0 and 0.0 are compared, yielding false.
			pathLength[i][j] = pathLength[i-1][j-1] + 1;
		else if (Double.compare(minGlobalCost + 0.0, fullICost + 0.0) == 0)
			pathLength[i][j] = pathLength[i-1][j] + 1;
		else if (Double.compare(minGlobalCost + 0.0, fullJCost + 0.0) == 0)
			pathLength[i][j] = pathLength[i][j-1] + 1;	

		normCostMatrix[i][j] = minGlobalCost;
	
	    } else {
		final double minGlobalCost = Math.min(costMatrix[i-1][j] + vcost + iCost,
		                                      Math.min(costMatrix[i-1][j-1] + diagCost,
		                                               costMatrix[i][j-1] + hcost + jCost));
		costMatrix[i][j] = minGlobalCost + distFn.calcDistance(tsI.getMeasurementVector(i),
		                                                           tsJ.getMeasurementVector(j));
	    }
	
         }  // end for loop
      }  // end for loop

//e = System.currentTimeMillis();
//System.out.println("Filling cost matrix: " + (e-s));

//s = System.currentTimeMillis();
      // Find the minimum score in the cost matrix in a cell (maxI, j). The entire query must be used, but the entire database need not be. 
      int bestJ = -1;
      double minimumCost = Double.POSITIVE_INFINITY;
      double[][] matrix = (normalize) ? normCostMatrix : costMatrix;
      for (int j = 0; j <= maxJ; j++) 
	if (matrix[maxI][j] < minimumCost)  {
		minimumCost = matrix[maxI][j];
		bestJ = j;
	}

//e = System.currentTimeMillis();
//System.out.println("Min: " + (e-s));

//s = System.currentTimeMillis();
      // Find the Warp Path by searching the matrix from the solution at
      //    (maxI, bestJ) up to the first cell (0,j), j <= bestJ.  At each step move through
      //    the matrix 1 step from (i, j) -> (i-1, j) or (i-1, j-1) or (i, j-1), whichever gives
      //    cell (i,j) the smallest cost.  Favour diagonal moves and moves towards the i==j
      //    axis to break ties.
      final WarpPath minCostPath = new WarpPath(maxI+maxJ-1);
      int i = maxI;
      int j = bestJ;
      minCostPath.addFirst(i, j);

      while (i>0)
      {
         // Find the costs of moving in all three possible directions.

	 final double diagCost; // motion weight cost for moving to i-1, j-1
         final double iCost; // motion weight cost for moving to i-1, j
         final double jCost; // motion weight cost for moving to i, j-1

         final double fullDiagCost; // full cost for moving to i-1, j-1
	 final double fullICost; // full cost for moving to i-1, j
	 final double fullJCost; // full cost for moving to i, j-1

	 double[] difI = tsIDiffs[i-1]; // The difference between frame i in tsI and frame i-1 in tsI.

 	 // Compute diagCost, iCost, and jCost, the motion weight costs for moving in the corresponding direction.
	 // Since we are in the loop, i > 0, so there will always be an iCost, and if j > 0, there will also be a diagCost and jCost. 
	 if (motionDist.equals("Euc")) {

		iCost = motionWeight*mot.calcDistance(tsIDiffs[i-1], zero);

		if (j > 0) {
			diagCost = motionWeight*mot.calcDistance(tsIDiffs[i-1], tsJDiffs[j-1]);
			jCost = motionWeight*mot.calcDistance(tsJDiffs[j-1], zero);
		} else {
			diagCost = 0; // It doesn't matter what these values are set as; they will not be used if j <= 0 (see "if (j > 0)" cases below).
			jCost = 0; 
		}

	} else if (motionDist.equals("Dot")) {

		// Compute iCost.
		double[] difJ;
		if (j == 0)
			difJ = tsJ.getMeasurementVector(0); // If j == 0, consider the frame before J's 0th frame to be the zero vector.
		else
			difJ = tsJDiffs[j-1];
		
		double[] jEpsilon = normalizeVector(momentumEpsilon(difJ));		
		iCost = motionWeight*mot.calcDistance(jEpsilon, difI);

		// Compute diagCost and jCost. 
		if (j > 0) {
			diagCost = motionWeight*mot.calcDistance(normalizeVector(tsIDiffs[i-1]), normalizeVector(tsJDiffs[j-1]));
			double[] iEpsilon = normalizeVector(momentumEpsilon(difI));		
			jCost = motionWeight*mot.calcDistance(difJ, iEpsilon);		
		} else {
			diagCost = 0; // It doesn't matter what these values are set as; they will not be used if j <= 0 (see "if (j > 0)" cases below).
			jCost = 0;
		}

	} else {
		diagCost = 0;
		iCost = 0;
		jCost = 0;
	}


	 // Compute the full costs for moving in each possible direction. 	
	 if (normalize) {
	         double myDist = distFn.calcDistance(tsI.getMeasurementVector(i), tsJ.getMeasurementVector(j));

		 fullICost = ( normCostMatrix[i-1][j]*pathLength[i-1][j] + iCost + vcost + myDist*motionWeight ) / //!!!
				(pathLength[i-1][j] + 1);
		 
		 if (j > 0) {
		    fullDiagCost = ( normCostMatrix[i-1][j-1]*pathLength[i-1][j-1] + diagCost + myDist ) /
				(pathLength[i-1][j-1] + 1);
		    fullJCost = ( normCostMatrix[i][j-1]*pathLength[i][j-1] + jCost + hcost + myDist*motionWeight) / //!!!
				(pathLength[i][j-1] + 1);
		 } else {//if j == 0, you can't move diagonally, or to (i, j-1)
		    fullDiagCost = Double.POSITIVE_INFINITY;
		    fullJCost = Double.POSITIVE_INFINITY;
		 }

	} else {
		 fullICost = costMatrix[i-1][j] + iCost + vcost;

		 if (j > 0) {
		    fullDiagCost = costMatrix[i-1][j-1] + diagCost;
		    fullJCost = costMatrix[i][j-1] + jCost + hcost;
		 } else {//if j == 0, you can't move diagonally, or to (i, j-1)
		    fullDiagCost = Double.POSITIVE_INFINITY;
		    fullJCost = Double.POSITIVE_INFINITY;
		 }
	}

         // Determine which direction to move in.  Prefer moving diagonally and
         //    moving towards the i==j axis of the matrix if there are ties.
         if ((fullDiagCost<=fullICost) && (fullDiagCost<=fullJCost))
         {
            i--;
            j--;
         }
         else if ((fullICost<fullDiagCost) && (fullICost<fullJCost))
            i--;
         else if ((fullJCost<fullDiagCost) && (fullJCost<fullICost))
            j--;
         else if (i <= j)  // iCost == jCost > diagCost
            j--;
         else   // iCost == jCost > diagCost
            i--;

         // Add the current step to the warp path.
         minCostPath.addFirst(i, j);
      }  // end while loop
//e = System.currentTimeMillis();
//System.out.println("Backtracking: " + (e-s));
//s = System.currentTimeMillis();
  //    TimeWarpInfo trimmed = trimCorners(new TimeWarpInfo(minimumCost, minCostPath), distFn, motionDist, motionWeight, vcost, hcost, tsIDiffs, tsJDiffs, tsI, tsJ, zero, 0);
//System.out.println("Backtracking: " + (e-s));
  //    return trimmed;
        return new TimeWarpInfo(minimumCost, minCostPath);   
}  // end DynamicTimeWarp(..)


   // These functions are the window-constrained versions of the ones above. The algorithm will search for the best warp path within a window given by the window parameter.
   public static WarpPath getWarpPathBetween(TimeSeries tsI, TimeSeries tsJ, SearchWindow window, DistanceFunction distFn)
   {
      return constrainedTimeWarp(tsI, tsJ, window, distFn, 0, 0, 0, "Dot", false).getPath();
   }

   public static WarpPath getWarpPathBetween(TimeSeries tsI, TimeSeries tsJ, SearchWindow window, DistanceFunction distFn, double vcost, double hcost, double motionWeight, String motionDist, boolean normalize)
   {
      return constrainedTimeWarp(tsI, tsJ, window, distFn, vcost, hcost, motionWeight, motionDist, normalize).getPath();
   }

   public static TimeWarpInfo getWarpInfoBetween(TimeSeries tsI, TimeSeries tsJ, SearchWindow window, DistanceFunction distFn)
   {
      return constrainedTimeWarp(tsI, tsJ, window, distFn, 0, 0, 0, "Dot", false);
   }

   public static TimeWarpInfo getWarpInfoBetween(TimeSeries tsI, TimeSeries tsJ, SearchWindow window, DistanceFunction distFn, double vcost, double hcost, double motionWeight, String motionDist, boolean normalize)
   {
      return constrainedTimeWarp(tsI, tsJ, window, distFn, vcost, hcost, motionWeight, motionDist, normalize);
   }

   // DTW, constrained to a given window of cells.
   private static TimeWarpInfo constrainedTimeWarp(TimeSeries tsI, TimeSeries tsJ, SearchWindow window, DistanceFunction distFn, double vcost, double hcost, double motionWeight, String motionDist, boolean normalize)
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

//double mt = 0;
//double ms;
//double me;
//double s;
//double e;
//s = System.currentTimeMillis();
      WindowMatrix costMatrix = null;
      WindowMatrix normCostMatrix = null;
      WindowMatrix pathLength = null;

      // For normalization, we need to keep track of normalized cost and path length, instead of total cost.
      if (normalize) {
	normCostMatrix = new WindowMatrix(window);
	pathLength = new WindowMatrix(window);
      } else {
	costMatrix = new WindowMatrix(window);
      }
//e = System.currentTimeMillis();
//System.out.println("Allocation: " + (e-s));

//s = System.currentTimeMillis();
      final int maxI = window.maxI();
      final int maxJ = window.maxJ();
      final int minJ = window.minJ();
      final int minI = 0; 
      final int dim = tsI.getMeasurementVector(0).length;
      // Create a zero vector, for calculating Euclidean distance motion weights when duplicating a frame (i,j) -> (i-1, j) or (i, j-1).
      final double[] zero = new double[dim];
      for (int z = 0; z < dim; z++)
         zero[z] = 0;

      // This is used for calculating distances for the motion weights.
      DistanceFunction mot = DistanceFunctionFactory.getDistFnByName(motionDist);
//e = System.currentTimeMillis();
//System.out.println("Misc: " + (e-s));

//s = System.currentTimeMillis();
      // Precompute the difference between each frame and the frame before; these values are used in the motion weights.  
      double[][] tsIDiffs = new double[maxI][dim];
      double[][] tsJDiffs = new double[maxJ - minJ + 1][dim];
      for (int tsind = 0; tsind < maxI; tsind++)
          tsIDiffs[tsind] = arrayMinus(tsI.getMeasurementVector(tsind+1), tsI.getMeasurementVector(tsind));
      for (int tsjind = minJ; tsjind <= maxJ; tsjind++) {
          if (tsjind != 0)
          	tsJDiffs[tsjind - minJ] = arrayMinus(tsJ.getMeasurementVector(tsjind), tsJ.getMeasurementVector(tsjind-1));
	  else
		tsJDiffs[tsjind - minJ] = tsJ.getMeasurementVector(tsjind); // If j == 0, consider the frame before J's 0th frame to be the zero vector.
      }
//e = System.currentTimeMillis();
//System.out.println("Velocities: " + (e-s));

      // Get an iterator that traverses the window cells in the order that the cost matrix is filled.
      //    (first to last row (1..maxI), bottom to top (1..MaxJ)
//s = System.currentTimeMillis();
      final Iterator matrixIterator = window.iterator();
//e = System.currentTimeMillis();
//System.out.println("Iterator: " + (e-s));

//s = System.currentTimeMillis();
      while (matrixIterator.hasNext())
      {
         final ColMajorCell currentCell = (ColMajorCell)matrixIterator.next();  // Current cell being filled.
         final int i = currentCell.getCol();
         final int j = currentCell.getRow();
	 
	 // If i == 0, no need to calculate any weights; just put the difference score between query frame i and database frame j into cell (0, j), since the warp path ends when all query frames are used.
         if (i == 0)
         {
	    if (normalize) {
		normCostMatrix.put(i, j, distFn.calcDistance(tsI.getMeasurementVector(0), tsJ.getMeasurementVector(j)));
		pathLength.put(i, j, 1);
	    } else
            	costMatrix.put(i, j, distFn.calcDistance(tsI.getMeasurementVector(0), tsJ.getMeasurementVector(j)));
            continue;
         }

	 /*motion weights costs for moving to [i-1][j], [i-1][j-1], or [j][j-1], respectively*/

//ms = System.currentTimeMillis();
	 double iCost = 0;
	 double diagCost = 0;
	 double jCost = 0;

	 // Here, i > 0 because of the conditional above, so set the motion weight cost of moving to (i-1, j), iCost.

	 switch (motionDist) {
		case "Euc":
		    iCost = motionWeight*mot.calcDistance(tsIDiffs[i-1], zero);
		    break;
	        case "Dot":
	            double[] difJ;
		    difJ = tsJDiffs[j-minJ];		
           	    double[] jEpsilon = normalizeVector(momentumEpsilon(difJ));		
		    iCost = motionWeight*mot.calcDistance(jEpsilon, normalizeVector(tsIDiffs[i-1]));
		    break;
	 }

	 // If j > 0, you can move to (i,j-1) or (i-1, j-1), since here, i > 0 because of the above conditional. 
	 if (j > 0) {
		switch (motionDist) {
			case "Euc":
				jCost = motionWeight*mot.calcDistance(tsJDiffs[j-minJ], zero);
				diagCost = motionWeight*mot.calcDistance(tsIDiffs[i-1], tsJDiffs[j-minJ]);
				break;
			case "Dot":
				double[] iEpsilon = normalizeVector(momentumEpsilon(tsIDiffs[i-1]));		
				jCost = motionWeight*mot.calcDistance(normalizeVector(tsJDiffs[j-minJ]), iEpsilon);	
				diagCost = 1*motionWeight*mot.calcDistance(normalizeVector(tsIDiffs[i-1]), normalizeVector(tsJDiffs[j-minJ]));	
				break;
		}
	 } else {
	 	jCost = 0; // If j == 0, jCost and and diagCost won't be used; it doesn't matter what they are set to.
	 	diagCost = 0;
	 }

/*	 if (motionDist.equals("Euc"))
		 iCost = motionWeight*mot.calcDistance(tsIDiffs[i-1], zero);
	 else if (motionDist.equals("Dot")) {
	         double[] difJ;
		 if (j == 0)
		 	difJ = tsJ.getMeasurementVector(0); //If j == 0, consider the frame before J's 0th frame to be the zero vector
		 else
			difJ = tsJDiffs[j-1];		
		 double[] jEpsilon = normalizeVector(momentumEpsilon(difJ));		
		 iCost = motionWeight*mot.calcDistance(jEpsilon, tsIDiffs[i-1]);
	 }

	 // If j > 0, you can move to (i,j-1) or (i-1, j-1), since here, i > 0 because of the above conditional. 
	 if (j > 0) {
		if (motionDist.equals("Euc")) {
			jCost = motionWeight*mot.calcDistance(tsJDiffs[j-1], zero);
			diagCost = motionWeight*mot.calcDistance(tsIDiffs[i-1], tsJDiffs[j-1]);
		} else if (motionDist.equals("Dot")) {
			double[] iEpsilon = normalizeVector(momentumEpsilon(tsIDiffs[i-1]));		
			jCost = motionWeight*mot.calcDistance(tsJDiffs[j-1], iEpsilon);	
			diagCost = 1*motionWeight*mot.calcDistance(normalizeVector(tsIDiffs[i-1]), normalizeVector(tsJDiffs[j-1]));	
		}
	 } else {
	 	jCost = 0; // If j == 0, jCost and and diagCost won't be used; it doesn't matter what they are set to.
	 	diagCost = 0;
	 }
*/
			
//me = System.currentTimeMillis();
//mt += (me-ms);       
        if (j == 0) // If j == 0, you can only move to (i-1, j).
         {
	    if (normalize) {
		 double fullICost = ( normCostMatrix.get(i-1,j)*pathLength.get(i-1,j) + motionWeight*distFn.calcDistance(tsI.getMeasurementVector(i), tsJ.getMeasurementVector(j)) + iCost + vcost ) / //!!!
					(pathLength.get(i-1,j) + 1); // This assumes that for window, if i1 < i2, then window.minJForI(i1) < window.minJforI(i2), since warp paths cannot move in directions of 									     // increasing j or i, and the window was calculated from a lower-resolution warp path.
		 normCostMatrix.put(i,j,fullICost);
		 pathLength.put(i, j, pathLength.get(i-1, j) + 1);
	    } else {
           	 costMatrix.put(i, j, distFn.calcDistance(tsI.getMeasurementVector(i), tsJ.getMeasurementVector(0)) +
                                 costMatrix.get(i-1, j) + iCost + vcost);
	    }

         }
         else // If j > 0, you can move to (i-1, j), (i-1, j-1), (i, j-1).
         {
	    if (normalize) {
		    double fullICost;
 		    double fullDiagCost;
 		    double fullJCost;
		    double pl;

		    double myDist = distFn.calcDistance(tsI.getMeasurementVector(i), tsJ.getMeasurementVector(j));

		    fullICost = ( normCostMatrix.get(i-1, j)*pathLength.get(i-1, j) + motionWeight*myDist + iCost + vcost ); //!!!
		    pl = pathLength.get(i-1,j) + 1;
		    if (pl != Double.POSITIVE_INFINITY) // If you are considering going onto a square outside of the window, don't divide by pl, or else your cost will become NaN.
			fullICost /= pl;

		    fullDiagCost = ( normCostMatrix.get(i-1,j-1)*pathLength.get(i-1, j-1) + myDist + diagCost);
		    pl = pathLength.get(i-1,j-1) + 1;
		    if (pl != Double.POSITIVE_INFINITY)
			fullDiagCost /= pl;

		    fullJCost = ( normCostMatrix.get(i,j-1)*pathLength.get(i,j-1) + motionWeight*myDist + jCost + hcost ); //!!!	   
		    pl = pathLength.get(i,j-1) + 1;
		    if (pl != Double.POSITIVE_INFINITY)
			fullJCost /= pl;

 		    final double minGlobalCost = Math.min(fullICost,
		                                          Math.min(fullJCost,
		                                                   fullDiagCost));

		    if (Double.compare(minGlobalCost + 0.0, fullICost + 0.0) == 0) {
			pathLength.put(i, j, pathLength.get(i-1,j) + 1);
		    } else if (Double.compare(minGlobalCost + 0.0, fullJCost + 0.0) == 0) {
		    	pathLength.put(i, j, pathLength.get(i,j-1) + 1);
		    } else if (Double.compare(minGlobalCost + 0.0, fullDiagCost + 0.0) == 0) {
			pathLength.put(i, j, pathLength.get(i-1,j-1) + 1);
		    }

		    normCostMatrix.put(i, j, minGlobalCost);
	    } else {
		    final double minGlobalCost = Math.min(costMatrix.get(i-1, j) + iCost + vcost,
		                                          Math.min(costMatrix.get(i-1, j-1) + diagCost,
		                                                   costMatrix.get(i, j-1) + jCost + hcost));
		    costMatrix.put(i, j, minGlobalCost + distFn.calcDistance(tsI.getMeasurementVector(i),
		                                                             tsJ.getMeasurementVector(j)));
	    }
         }  // end if
      }  // end while loop
//e = System.currentTimeMillis();
//System.out.println("Filling CM: " + (e-s) + "  mt: " + mt);

//s = System.currentTimeMillis();
      // Find the minimum score in the cost matrix in a cell (maxI, j). The entire query must be used, but the entire database need not be. 
      int bestJ = -1;
      double minimumCost = Double.POSITIVE_INFINITY;
      WindowMatrix matrix = (normalize) ? normCostMatrix : costMatrix;
      for (int j = window.minJforI(maxI); j <= maxJ; j++) 
	if (matrix.get(maxI, j) < minimumCost)  {
		minimumCost = matrix.get(maxI, j);
		bestJ = j;
	}
//e = System.currentTimeMillis();
//System.out.println("Finding bestJ: " + (e-s));

//s = System.currentTimeMillis();
      // Find the Warp Path by searching the matrix from the solution at
      //    (maxI, bestJ) up to the first cell (0,j), j <= bestJ.  At each step move through
      //    the matrix 1 step from (i, j) -> (i-1, j) or (i-1, j-1) or (i, j-1), whichever gives
      //    cell (i,j) the smallest cost.  Favour diagonal moves and moves towards the i==j
      //    axis to break ties.
      final WarpPath minCostPath = new WarpPath(maxI+maxJ-minJ); //maxI + maxJ - minJ is the longest warp path possible, obtained by not moving diagonally
      int i = maxI;
      int j = bestJ;
      minCostPath.addFirst(i, j);
      while (i>0)
      {
         // Find the costs of moving in all three possible directions.
	 
         final double diagCost; // motion weight cost for moving to i-1, j-1
         final double iCost; // motion weight cost for moving to i-1, j
         final double jCost; // motion weight cost for moving to i, j-1

	 double fullDiagCost; // full cost for moving to i-1, j-1
	 double fullICost; // full cost for moving to i-1, j
	 double fullJCost; // full cost for moving to i, j-1

	 final double[] difI = tsIDiffs[i-1]; // The difference between frame i in tsI and frame i-1 in tsI.
	 double[] difJ;

	 // Compute diagCost, iCost, and jCost, the motion weight costs for moving in the corresponding direction.
	 // Since we are in the loop, i > 0, so there will always be an iCost, and if j > 0, there will also be a diagCost and jCost. 

	 if (motionDist.equals("Euc"))
      	        iCost = motionWeight*mot.calcDistance(difI, zero);       
	 else if (motionDist.equals("Dot")) {
		difJ = tsJDiffs[j-minJ];		
		double[] jEpsilon = normalizeVector(momentumEpsilon(difJ));		
		iCost = motionWeight*mot.calcDistance(jEpsilon, difI);	  
	 } else
	        iCost = 0;

	 // Compute the motion weight cost for moving to i, j - 1, and for moving to i - 1, j -1
         if (j > 0) {
	    difJ = tsJDiffs[j-minJ];

	    if (motionDist.equals("Euc"))
            	diagCost = motionWeight*mot.calcDistance(difI, difJ);
	    else if (motionDist.equals("Dot"))
		diagCost = 1*motionWeight*mot.calcDistance(normalizeVector(difI), normalizeVector(difJ));
	    else diagCost = 0;

	    if (motionDist.equals("Euc"))
	    	jCost = motionWeight*mot.calcDistance(difJ, zero);
	    else if (motionDist.equals("Dot")) {
		double[] iEpsilon = normalizeVector(momentumEpsilon(tsIDiffs[i-1]));		
		jCost = motionWeight*mot.calcDistance(tsJDiffs[j-minJ], iEpsilon);
	    } else
		jCost = 0;

	 } else {
            diagCost = Double.POSITIVE_INFINITY;
	    jCost = Double.POSITIVE_INFINITY; 
	 }


	 if (normalize) {
		double myDist = distFn.calcDistance(tsI.getMeasurementVector(i), tsJ.getMeasurementVector(j));
                double pl; //pathLength
		//Compute the full cost for moving to i-1,j
		fullICost = ( normCostMatrix.get(i-1, j)*pathLength.get(i-1, j) + iCost + vcost + motionWeight*myDist );  //!!!
		pl = (pathLength.get(i-1, j) + 1);
		if (pl != Double.POSITIVE_INFINITY)
			fullICost = fullICost / pl;

		//Compute the full cost for moving to i-1,j-1 and i, j-1
		if (j > 0) {
			fullJCost = ( normCostMatrix.get(i, j-1)*pathLength.get(i, j-1) + jCost + hcost + motionWeight*myDist ); //!!!
			pl = (pathLength.get(i, j-1) + 1);
			if (pl != Double.POSITIVE_INFINITY)
				fullJCost = fullJCost / pl;
			fullDiagCost = ( normCostMatrix.get(i-1, j-1)*pathLength.get(i-1, j-1) + diagCost + myDist );
			pl = (pathLength.get(i-1, j-1) + 1);
			if (pl != Double.POSITIVE_INFINITY)
				fullDiagCost = fullDiagCost / pl;
		} else {
		  	fullDiagCost = Double.POSITIVE_INFINITY;
			fullJCost = Double.POSITIVE_INFINITY;
		}

	 } else {
		//Compute the full cost for moving to i-1,j
		fullICost = costMatrix.get(i-1, j) + vcost + iCost;

		//Compute the full cost for moving to i-1,j-1 and i, j-1
		if (j > 0) {
			fullDiagCost = costMatrix.get(i-1, j-1) + diagCost;
			fullJCost = costMatrix.get(i, j-1) + hcost + jCost;
		} else {
			fullDiagCost = Double.POSITIVE_INFINITY;
			fullJCost = Double.POSITIVE_INFINITY;
		}
	 }

         // Determine which direction to move in.  Prefer moving diagonally and
         //    moving towards the i==j axis of the matrix if there are ties.
         if ((fullDiagCost<=fullICost) && (fullDiagCost<=fullJCost))
         {
            i--;
            j--;
         }
         else if ((fullICost<fullDiagCost) && (fullICost<fullJCost))
            i--;
         else if ((fullJCost<fullDiagCost) && (fullJCost<fullICost))
            j--;
         else if (i <= j)  // leftCost==rightCost > diagCost
            j--;
         else   // leftCost==rightCost > diagCost
            i--;

         // Add the current step to the warp path.
         minCostPath.addFirst(i, j);
      }  // end while loop
//e = System.currentTimeMillis();
//System.out.println("Backtracking: " + (e-s));

//s = System.currentTimeMillis();
      // Free any rescources associated with the costMatrix (a swap file may have been created if the swa file did not
      //    fit into main memory).
      if (costMatrix != null)
    	  costMatrix.freeMem();
//e = System.currentTimeMillis();
//System.out.println("Freeing memory: " + (e-s));

//s = System.currentTimeMillis();
    //  TimeWarpInfo trimmed = trimCorners(new TimeWarpInfo(minimumCost, minCostPath), distFn, motionDist, motionWeight, vcost, hcost, tsIDiffs, tsJDiffs, tsI, tsJ, zero, minJ);
//e = System.currentTimeMillis();
//System.out.println("Trimming: " + (e-s));
 //     return trimmed;

	return new TimeWarpInfo(minimumCost, minCostPath);
   }  // end ConstrainedTimeWarp


   // Helper for substracting two vectors. 
   private static double[] arrayMinus(double[] b, double[] a) {
	if (a.length != b.length)
		return null; 
	double [] c = new double[a.length];	
	for (int i = 0; i < a.length; i++)
		c[i] = b[i] - a[i];
        return c;
   }

//public static double normTime = 0;
   // Helper for normalizing a vector.
   private static double[] normalizeVector(double[] v) {
//double s = System.currentTimeMillis();
	double norm = 0;	
	for (int i = 0; i < v.length; i++) {
		norm += v[i]*v[i];
	}

	norm = Math.sqrt(norm);
        
        double[] r = new double[v.length];

        if (norm == 0) {
           for (int i = 0; i < v.length; i++) {
              r[i] = 0;
           }
        } else {
	   for (int i = 0; i < v.length; i++) {
	      r[i] = v[i] / norm;
	   }
        }
//double e = System.currentTimeMillis();
//normTime += (e-s);
	return r;
   }

   // Helper for calculating the "epsilon" from a vector, which is the diagonal vector in the quadrant that the input vector is part of. The idea is that this vector should
   // have a higher dot product with a vector representing small frame to frame motion, and a lower dot product with a vector representing large motion.
   private static double[] momentumEpsilon(double[] v) {
	double[] r = new double[v.length];
	for (int i = 0; i < v.length; i++) 
		if (v[i] < 0)
			r[i] = -1;
		else
			r[i] = 1;
	return r;
   }

   // Helper for trimming a warp path of corners _|, because these seem weird when interpreting the match.
   // zero is the zero vector of dimensionality of the time series
   // dFn is the Distance function for frame-to-frame comparison. motionDist specifies
   // the Distance function for motion weights (ie. frame-velocity to frame-velocity comparison)
   private static TimeWarpInfo trimCorners(TimeWarpInfo twp, DistanceFunction dFn, String motionDist, double motionWeight, double vcost, double hcost, double[][] tsIDiffs, double[][] tsJDiffs, TimeSeries tsI, TimeSeries tsJ, double[] zero, int minJ) {
       WarpPath path = twp.getPath();
       double score = twp.getDistance()*path.size();

       // You can't have corners if your path is puny
       if (path.size() < 3) 	
	  return twp;

       DistanceFunction mot = DistanceFunctionFactory.getDistFnByName(motionDist);
       
       // initialize
       int lasti2 = path.get(0).getCol();
       int lastj2 = path.get(0).getRow();
       int lasti = path.get(1).getCol();
       int lastj = path.get(1).getRow();

       int w = 2;
       while (w < path.size()) {
           ColMajorCell c = path.get(w);
           int i = c.getCol();
           int j = c.getRow();
           
           int idif = i - lasti;
           int jdif = j - lastj;
	   int idif2 = lasti - lasti2;
	   int jdif2 = lastj - lastj2;

           // "If you encounter a corner, trim it and update the score..."
           if ((idif == 0 && jdif == 1 && idif2 == 1 && jdif2 == 0) ||
               (idif == 1 && jdif == 0 && idif2 == 0 && jdif2 == 1)) {
	        
	   	// Add the motion weight for diagonal movement.
                double[] diagIVelocity = tsIDiffs[i-1];
		double[] diagJVelocity = tsJDiffs[j-minJ];

                if (motionDist.equals("Euc")) {
			score += motionWeight*mot.calcDistance(diagIVelocity, diagJVelocity); // Don't normalize with Euclidean distance motion costs, as this
									        	      // is intended to compare quantity of movement as well 
                } else if (motionDist.equals("Dot")) {
			score += motionWeight*mot.calcDistance(normalizeVector(diagIVelocity), normalizeVector(diagJVelocity));
                }
   
                // Subtract the motion weights for the vertical and horizontnal movements.
		score -= vcost;
		score -= hcost;
		if (motionDist.equals("Euc")) {
			score -= motionWeight*mot.calcDistance(diagIVelocity, zero);
			score -= motionWeight*mot.calcDistance(diagJVelocity, zero);
		} else if (motionDist.equals("Dot")) {
		        if (idif == 1) { // corner: first move in j direction, then in i direction
				double[] iEpsilon = momentumEpsilon(tsIDiffs[i-2]); //This contains the velocity moving from frame i-2 to i-1
				score -= motionWeight*mot.calcDistance(normalizeVector(iEpsilon), normalizeVector(diagJVelocity));
				double[] jEpsilon = momentumEpsilon(diagJVelocity); 
				score -= motionWeight*mot.calcDistance(normalizeVector(jEpsilon), normalizeVector(diagIVelocity));

		        } else { // corner: first move in i direction, then in j direction
				double[] iEpsilon = momentumEpsilon(diagIVelocity);
				score -= motionWeight*mot.calcDistance(normalizeVector(iEpsilon), normalizeVector(diagJVelocity));
				double[] jEpsilon = momentumEpsilon(tsJDiffs[j-minJ-1]); //This contains the velocity moving from frame j-2 to j-1
				score -= motionWeight*mot.calcDistance(normalizeVector(jEpsilon), normalizeVector(diagIVelocity));
		        }
		}

		// Subtract the score for the additional frame.
		score -= motionWeight*dFn.calcDistance(tsI.getMeasurementVector(lasti), tsJ.getMeasurementVector(lastj));

		// Remove the corner from the warp path
		path.remove(w-1);

		// Update the indices
		lasti = i;
		lastj = j;

 	   } else {
	        // Increment w if you're not removing anything. If you are, no need, as the next item will shift down by 1.
		w++;

		lasti2 = lasti;
		lastj2 = lastj;
		lasti = i;
		lastj = j;
	   }
       }

    score /= path.size();

    return new TimeWarpInfo(score, path);
	
   } 

}  // end class DtwTest

/*
 * FastDTW.java   Jul 14, 2004
 *
 * Copyright (c) 2004 Stan Salvador
 * stansalvador@hotmail.com
 */

package com.dtw;

import com.timeseries.TimeSeries;
import com.timeseries.PAA;
import com.util.DistanceFunction;

import java.lang.String;


public class FastDTW
{
   // CONSTANTS
   final static int DEFAULT_SEARCH_RADIUS = 1;

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

   /* Basic functions to call Fast DTW, setting only the main options (Distance Function ("DotProductDistance", etc.), and the window radius for FastDTW) */
   public static double getWarpDistBetween(TimeSeries tsI, TimeSeries tsJ, DistanceFunction distFn)
   {
      return fastDTW(tsI, tsJ, DEFAULT_SEARCH_RADIUS, distFn, 0, 0, 0, "Dot", false).getDistance();
   }


   public static double getWarpDistBetween(TimeSeries tsI, TimeSeries tsJ, int searchRadius, DistanceFunction distFn)
   {
      return fastDTW(tsI, tsJ, searchRadius, distFn, 0, 0, 0, "Dot", false).getDistance();
   }


   public static WarpPath getWarpPathBetween(TimeSeries tsI, TimeSeries tsJ, DistanceFunction distFn)
   {
      return fastDTW(tsI, tsJ, DEFAULT_SEARCH_RADIUS, distFn, 0, 0, 0, "Dot", false).getPath();
   }


   public static WarpPath getWarpPathBetween(TimeSeries tsI, TimeSeries tsJ, int searchRadius, DistanceFunction distFn)
   {
      return fastDTW(tsI, tsJ, searchRadius, distFn, 0, 0, 0, "Dot", false).getPath();
   }

 
   public static TimeWarpInfo getWarpInfoBetween(TimeSeries tsI, TimeSeries tsJ, int searchRadius, DistanceFunction distFn)
   {
      return fastDTW(tsI, tsJ, searchRadius, distFn, 0, 0, 0, "Dot", false);
   }


  /* Versions of the above functions with additional parameters */
  public static double getWarpDistBetween(TimeSeries tsI, TimeSeries tsJ, DistanceFunction distFn, double vcost, double hcost, double motionWeight, String motionDist, boolean normalize)
   {
      return fastDTW(tsI, tsJ, DEFAULT_SEARCH_RADIUS, distFn, vcost, hcost, motionWeight, motionDist, normalize).getDistance();
   }


   public static double getWarpDistBetween(TimeSeries tsI, TimeSeries tsJ, int searchRadius, DistanceFunction distFn, double vcost, double hcost, double motionWeight, String motionDist, boolean normalize)
   {
      return fastDTW(tsI, tsJ, searchRadius, distFn, vcost, hcost, motionWeight, motionDist, normalize).getDistance();
   }


   public static WarpPath getWarpPathBetween(TimeSeries tsI, TimeSeries tsJ, DistanceFunction distFn, double vcost, double hcost, double motionWeight, String motionDist, boolean normalize)
   {
      return fastDTW(tsI, tsJ, DEFAULT_SEARCH_RADIUS, distFn, vcost, hcost, motionWeight, motionDist, normalize).getPath();
   }


   public static WarpPath getWarpPathBetween(TimeSeries tsI, TimeSeries tsJ, int searchRadius, DistanceFunction distFn, double vcost, double hcost, double motionWeight, String motionDist, boolean normalize)
   {
      return fastDTW(tsI, tsJ, searchRadius, distFn, vcost, hcost, motionWeight, motionDist, normalize).getPath();
   }


   public static TimeWarpInfo getWarpInfoBetween(TimeSeries tsI, TimeSeries tsJ, int searchRadius, DistanceFunction distFn, double vcost, double hcost, double motionWeight, String motionDist, boolean normalize)
   {
      return fastDTW(tsI, tsJ, searchRadius, distFn, vcost, hcost, motionWeight, motionDist, normalize);
   }


   /* The recursive Fast DTW Function */
   private static TimeWarpInfo fastDTW(TimeSeries tsI, TimeSeries tsJ, int searchRadius, DistanceFunction distFn, double vcost, double hcost, double motionWeight, String motionDist, boolean normalize)
   {
      if (searchRadius < 0)
         searchRadius = 0;

      final int minTSsize = searchRadius+2;

      // Base Case
      if ( (tsI.size() <= minTSsize) || (tsJ.size()<=minTSsize) )
      {
         // Perform full Dynamic Time Warping.
         TimeWarpInfo info = DTW.getWarpInfoBetween(tsI, tsJ, distFn, vcost, hcost, motionWeight, motionDist, normalize);
         return info;
      }
      else
      {
         final double resolutionFactor = 2.0;

         final PAA shrunkI = new PAA(tsI, (int)(tsI.size()/resolutionFactor));
         final PAA shrunkJ = new PAA(tsJ, (int)(tsJ.size()/resolutionFactor));

 
        // Determine the search window that constrains the area of the cost matrix that will be evaluated based on
        //    the warp path found at the previous resolution (smaller time series).
	WarpPath shrunkWarpPath = FastDTW.getWarpPathBetween(shrunkI, shrunkJ, searchRadius, distFn, vcost, hcost, motionWeight, motionDist, normalize);

         final SearchWindow window = new ExpandedResWindow(tsI, tsJ, shrunkI, shrunkJ,
                                                            shrunkWarpPath,
                                                            searchRadius);

         // Find the optimal warp path through this search window constraint.
         return DTW.getWarpInfoBetween(tsI, tsJ, window, distFn, vcost, hcost, motionWeight, motionDist, normalize);
      }  // end if
   }  // end recFastDTW(...)

}  // end class fastDTW

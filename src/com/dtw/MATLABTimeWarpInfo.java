
package com.dtw;

import com.dtw.TimeWarpInfo;

public class MATLABTimeWarpInfo
{
   // PRIVATE DATA
   private final double distance;
   private final double[][] path;


   // CONSTRUCTOR
   public MATLABTimeWarpInfo(TimeWarpInfo info)
   {
      distance = info.getDistance();
      path = info.getPath().toMATLAB();
   }

   public double getDistance() {
	return distance;
   }

   public double[][] getPath() {
	return path;
   }

}  // end class MATLABTimeWarpInfo

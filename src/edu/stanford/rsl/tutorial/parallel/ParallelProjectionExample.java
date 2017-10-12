package edu.stanford.rsl.tutorial.parallel;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.Grid2DComplex;
import edu.stanford.rsl.tutorial.filters.RamLakKernel;
import edu.stanford.rsl.tutorial.phantoms.DotsGrid2D;
import edu.stanford.rsl.tutorial.phantoms.MickeyMouseGrid2D;
import edu.stanford.rsl.tutorial.phantoms.Phantom;
import edu.stanford.rsl.tutorial.phantoms.TestObject1;
import edu.stanford.rsl.tutorial.phantoms.UniformCircleGrid2D;
import ij.ImageJ;

/**
 * This example can be used to demonstrate projections of different phantoms.
 * 
 * @author Andreas Maier
 *
 */
public class ParallelProjectionExample {

	public static void main (String [] args){
		new ImageJ();
		
		int x = 512;
		int y = 512;
		// Create a phantom
//		Phantom phan = new TestObject1(x, y);
//		phan = new UniformCircleGrid2D(x, y);
//		phan = new MickeyMouseGrid2D(x, y, 3);
		
		Phantom phan = new TestObject1(512, 512);
		phan.show("The Phantom");
		
		// Project forward parallel
		ParallelProjector2D projector = new ParallelProjector2D(Math.PI, Math.PI/1200.0, 1000, 0.5);
		Grid2D sinogram = projector.projectRayDrivenCL(phan);
		sinogram.show("The Sinogram");
		
		
		Grid2DComplex sino_fft = new Grid2DComplex(sinogram);
		sino_fft.transformForward();
		sino_fft.fftshift();
		sino_fft.show("fft of sinogram");
		
		
		
	}
	
}
/*
 * Copyright (C) 2014 Andreas Maier
 * CONRAD is developed as an Open Source project under the GNU General Public License (GPL).
*/
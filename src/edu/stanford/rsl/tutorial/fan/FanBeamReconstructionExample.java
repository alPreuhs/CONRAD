package edu.stanford.rsl.tutorial.fan;

import ij.ImageJ;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.Grid2DComplex;
import edu.stanford.rsl.conrad.data.numeric.NumericGrid;
import edu.stanford.rsl.conrad.data.numeric.NumericPointwiseOperators;
import edu.stanford.rsl.conrad.filtering.redundancy.ParkerWeightingTool;
import edu.stanford.rsl.conrad.utils.CONRAD;
import edu.stanford.rsl.tutorial.fan.FanBeamBackprojector2D;
import edu.stanford.rsl.tutorial.fan.FanBeamProjector2D;
import edu.stanford.rsl.tutorial.fan.redundancy.BinaryWeights;
import edu.stanford.rsl.tutorial.fan.redundancy.CompensationWeights;
import edu.stanford.rsl.tutorial.fan.redundancy.ParkerWeights;
import edu.stanford.rsl.tutorial.fan.redundancy.SilverWeights;
import edu.stanford.rsl.tutorial.filters.RamLakKernel;
import edu.stanford.rsl.tutorial.phantoms.DotsGrid2D;
import edu.stanford.rsl.tutorial.phantoms.MickeyMouseGrid2D;
import edu.stanford.rsl.tutorial.phantoms.Phantom;
import edu.stanford.rsl.tutorial.phantoms.TestObject1;
import edu.stanford.rsl.tutorial.phantoms.UniformCircleGrid2D;

/**
 * Simple example that computes and displays a reconstruction.
 * 
 * @author Recopra Seminar Summer 2012
 * 
 */
public class FanBeamReconstructionExample {


	/**
	 * @param args
	 */
	public static void main(String[] args) {
		CONRAD.setup();
		
		ParkerWeightingTool pwt = new ParkerWeightingTool();
	
		try {
			pwt.configure();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// image params
		int imgSzXMM = 512, // [mm]
		imgSzYMM = imgSzXMM; // [mm]
		float pxSzXMM = 1.0f, // [mm]
		pxSzYMM = pxSzXMM; // [mm]
		// fan beam bp parameters
		
		int numProj = 1600;
		double focalLength = 1350;
		double maxT =700;
		double deltaT = 1.0;
		double gammaM = Math.atan((maxT/2.0)/focalLength);
		double maxBeta = Math.PI+2*gammaM;
		double deltaBeta = maxBeta/numProj;
		

		///// compute primary angle
		double[] primary_angles = new double[numProj];
		for(int i =0;i <primary_angles.length; i++) {
			primary_angles[i] = Math.toDegrees(deltaBeta*i);				
		}

		pwt.setNumberOfProjections(numProj);
		pwt.setDetectorWidth((int)maxT);
		pwt.setPixelDimensionX(deltaT);
		pwt.setSourceToDetectorDistance(focalLength);
		pwt.setPrimaryAngles(primary_angles);
		
		
		Grid2D pw = new Grid2D((int) maxT,numProj);
		
		
		for(int i = 0 ; i< numProj; i++) {
			double[] parker = pwt.computeParkerWeights1D(i);
			for(int j = 0 ; j<parker.length;j++) {
				pw.getSubGrid(i).setAtIndex(j, (float) parker[j]); ;
			}
		}
		pw.show("parker weights using pwt");
		

		int phantomType = 2; // 0 = circle, 1 = MickeyMouse, 2 = TestObject1,
		// 3=DotsGrid
		// size in grid unit
		
		int imgSzXGU = (int) Math.floor(imgSzXMM / pxSzXMM), // [GU]
		imgSzYGU = (int) Math.floor(imgSzYMM / pxSzYMM); // [GU]
		new ImageJ();

		FanBeamProjector2D fanBeamProjector = new FanBeamProjector2D(
				focalLength, maxBeta, deltaBeta, maxT, deltaT);

		Phantom phantom;
		switch (phantomType) {
		case 0:
			phantom = new UniformCircleGrid2D(imgSzXGU, imgSzYGU);
			break;
		case 1:
			phantom = new MickeyMouseGrid2D(imgSzXGU, imgSzYGU);
			break;
		case 2:
			phantom = new TestObject1(imgSzXGU, imgSzYGU);
			break;
		case 3:
			phantom = new DotsGrid2D(imgSzXGU, imgSzYGU);
			break;
		default:
			phantom = new UniformCircleGrid2D(imgSzXGU, imgSzYGU);
			break;
		}

		phantom.setSpacing(pxSzXMM, pxSzYMM);
		// origin is given in (negative) world coordinates
		phantom.setOrigin(-(imgSzXGU * phantom.getSpacing()[0]) / 2, -(imgSzYGU * phantom.getSpacing()[1]) / 2);
		//phantom.setOrigin(-50.0, -50.0);

		phantom.show();
		Grid2D projectionP = new Grid2D(phantom);		
		
		for (int iter =0; iter < 1; iter ++) {
			// create projections
			Grid2D fanBeamSinoRay = fanBeamProjector.projectRayDrivenCL(projectionP);
			
			Grid2DComplex sino_fft = new Grid2DComplex(fanBeamSinoRay);
			sino_fft.transformForward();
			sino_fft.fftshift();
			sino_fft.show("fft of sinogram");
			
			
			fanBeamSinoRay.clone().show("Sinogram");

			int weightType = 2;
			
			Grid2D RedundancyWeights;
			switch (weightType) {
			case 0:
				RedundancyWeights = new ParkerWeights(focalLength, maxT, deltaT, maxBeta, deltaBeta);
				break;
			case 1:
				RedundancyWeights = new SilverWeights(focalLength, maxT, deltaT, maxBeta, deltaBeta);
				break;
			case 2:
				RedundancyWeights = new CompensationWeights(focalLength, maxT, deltaT, maxBeta, deltaBeta);
				break;
			case 3:
				RedundancyWeights = new BinaryWeights(focalLength, maxT, deltaT, maxBeta, deltaBeta);
				break;
			default:
				RedundancyWeights = new CompensationWeights(focalLength, maxT, deltaT, maxBeta, deltaBeta);
			}
			
			RedundancyWeights.show("Current Weight");

			RamLakKernel ramLak = new RamLakKernel((int) (maxT / deltaT), deltaT);
			CosineFilter cKern = new CosineFilter(focalLength, maxT, deltaT);
			// Apply filtering
			for (int theta = 0; theta < fanBeamSinoRay.getSize()[1]; ++theta) {
				cKern.applyToGrid(fanBeamSinoRay.getSubGrid(theta));
				ramLak.applyToGrid(fanBeamSinoRay.getSubGrid(theta));
			}
			NumericPointwiseOperators.multiplyBy(fanBeamSinoRay, pw);
			
			fanBeamSinoRay.show("After Filtering");
			
			// Do the backprojection
			FanBeamBackprojector2D fbp = new FanBeamBackprojector2D(focalLength,
					deltaT, deltaBeta, imgSzXMM, imgSzYMM);

			Grid2D reco = fbp.backprojectPixelDrivenCL(fanBeamSinoRay);
			reco.show("Parker" + iter);

			NumericGrid recoDiff = NumericPointwiseOperators.subtractedBy(phantom, reco);
			recoDiff.show("RecoDiff" + iter);
		}
	}
}
/*
 * Copyright (C) 2010-2014 Andreas Maier
 * CONRAD is developed as an Open Source project under the GNU General Public License (GPL).
*/
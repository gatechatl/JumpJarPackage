package PurityEstimation;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.LinkedList;

import org.systemsbiology.jrap.DataProcessingInfo;
import org.systemsbiology.jrap.MSInstrumentInfo;
import org.systemsbiology.jrap.MSXMLParser;
import org.systemsbiology.jrap.MZXMLFileInfo;
import org.systemsbiology.jrap.Scan;
import org.openscience.cdk.formula.IsotopeContainer;
import org.openscience.cdk.formula.IsotopePattern;
import org.openscience.cdk.formula.IsotopePatternManipulator;

import CDKFunction.IsotopePatternSimilarity;
import METABOLOMICS.AIM.AIMPuritySampler;

public class IsotopePurityMeasurement {
	private static double C12 = 0.9893;
	private static double N14 = 0.99632;	
	private static double H1 = 0.999885;		
	public static String description() {
		return "Parse a particular scan and generate the peak list";
	}
	public static String type() {
		return "MZXML";
	}
	public static String parameter_info() {
		return "[inputmzXMLFile] [scanNumber] [formula] [monoisotopicmass] [charge] [ppm]";
	}
	public static void main(String[] args) {
		System.out.println("mzXMLParser");
	}
	public static void execute(String[] args) {
		
		try {
			
			String inputFile = args[0];			
			int scanNumber = new Integer(args[1]);
			String formula = args[2];
			double monoIsotopicMass = new Double(args[3]);
			int charge = new Integer(args[4]);
			double ppm = new Double(args[5]);
			
			MSXMLParser parser = new MSXMLParser(inputFile);
			Scan scan = parser.rap(scanNumber);
			float[][] peaks = scan.getMassIntensityList();
			MZXMLFileInfo info = parser.getHeaderInfo();
			DataProcessingInfo dataProcessingInfo = info.getDataProcessing();
			IsotopePattern query_peaks = MatchIsotopeToolbox.getPeakInfo(peaks);
			
			IsotopePattern clean_query_peaks = MatchIsotopeToolbox.cleanPeak(query_peaks, 0.02); // clean the peaks

			for (int i = 0; i < clean_query_peaks.getNumberOfIsotopes(); i++) {
				IsotopeContainer container = clean_query_peaks.getIsotope(i);
				double mass = container.getMass();
				double intensity = container.getIntensity();
				//System.out.println(mass + "\t" + intensity);
			}
			
			MSInstrumentInfo instrumentInfo = info.getInstrumentInfo();
			//System.out.println("Instrument Info: " + instrumentInfo.getModel());
			//System.out.println(peaks[0].length);
			//System.out.println(peaks[1].length);
			for (int i = 0; i < clean_query_peaks.getNumberOfIsotopes(); i++) {
				IsotopeContainer container = clean_query_peaks.getIsotope(i);
				double mass = container.getMass();
				double intensity = container.getIntensity();
				//System.out.println("cleanedPeaks\t" + mass + "\t" + intensity);
			}
			//System.out.println("### End reading peaks ###");
			double minMass = Double.MAX_VALUE;
			double maxMass = Double.MIN_VALUE;
			AIMPuritySampler sampler = new AIMPuritySampler();
			IsotopePattern theoretical_isotope = sampler.calculate_pattern(formula, C12, N14, H1, charge);
			
			//System.out.println("### End Generating Theoretical Isotope ###");
			
			IsotopePattern clean_query_peaks_target = MatchIsotopeToolbox.minimizeArea(clean_query_peaks, theoretical_isotope, ppm);
			IsotopePattern clean_query_peaks_target_norm = IsotopePatternManipulator.normalize(clean_query_peaks_target);
			
			IsotopePattern c12_removed_clean_query_peaks_target = MatchIsotopeToolbox.removeC12Peak(clean_query_peaks, clean_query_peaks_target);
			
			//IsotopeContainer container = theoretical_isotope.getIsotope(0);
			
			
			IsotopePatternSimilarity similarity = new IsotopePatternSimilarity();
			double similarity_score = similarity.compare(theoretical_isotope, clean_query_peaks_target_norm);
			System.out.println("Similarity_score: " + similarity_score);
			
			// find C13 peaks
			double[] c12_simulation_result = MatchIsotopeToolbox.simulate_c13(formula, clean_query_peaks_target_norm, ppm, charge, 0, 1);			
			double[] c13_simulation_result = MatchIsotopeToolbox.simulate_c13(formula, c12_removed_clean_query_peaks_target, ppm, charge, 0, 1);			
			double[] n15_simulation_result = MatchIsotopeToolbox.simulate_n15(formula, c12_removed_clean_query_peaks_target, ppm, charge, 0, 1);
			
			System.out.println("Final Hit: " + c13_simulation_result[0] + "\t" + c13_simulation_result[1]);
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	/*public static double[] simulate_n15(String formula, IsotopePattern c12_removed_clean_query_peaks_target, double ppm, int charge) {
		AIMPuritySampler sampler = new AIMPuritySampler();
		// find C13 peaks
		double big_step = 0.05;
		boolean doOnce = true;
		double smallest = 0;
		double largest = 1;
		double prevIndex = 0;
		double step = 0.01;
		double simulated_N15_max_score = -1;
		IsotopePatternSimilarity similarity = new IsotopePatternSimilarity();
		for (double N14_sample = 0; N14_sample <= 1; N14_sample += big_step) {
			IsotopePattern theoretical_isotope_n14 = sampler.calculate_pattern(formula, C12, N14_sample, H1, charge);
			IsotopePattern test_target = minimizeArea(c12_removed_clean_query_peaks_target, theoretical_isotope_n14, ppm);
			IsotopePattern test_target_norm = IsotopePatternManipulator.normalize(test_target);
			if (test_target_norm.getNumberOfIsotopes() > 0) {
				double score = similarity.compare(test_target_norm, theoretical_isotope_n14);
				if (score > 0) {
					if (doOnce) {
						smallest = prevIndex;
						doOnce = false;
					}
					largest = N14_sample;
				}
			}
			prevIndex = N14_sample;
		}
		largest = largest + big_step;
		if (largest <= 1 - big_step) {
			largest += big_step;
			if (largest >= 1) {
				largest = 1.0;
			}
		}
		double max_score = 0;
		for (double N14_sample = 0; N14_sample < smallest; N14_sample += step) {
			//System.out.println(C12_sample + "\t" + 0.0);
		}
		for (double N14_sample = smallest; N14_sample <= largest; N14_sample += step) {			
			IsotopePattern theoretical_isotope_c13 = sampler.calculate_pattern(formula, C12, N14_sample, H1, charge);
			IsotopePattern test_target = minimizeArea(c12_removed_clean_query_peaks_target, theoretical_isotope_c13, ppm);
			IsotopePattern test_target_norm = IsotopePatternManipulator.normalize(test_target);
			if (test_target_norm.getNumberOfIsotopes() > 0) {
				double score = similarity.compare(test_target_norm, theoretical_isotope_c13);
				if (max_score < score) {
					max_score = score;
					simulated_N15_max_score = N14_sample;
				}
				
				//System.out.println(C12_sample + "\t" + score);
			}
		}
		for (double N14_sample = largest + step; N14_sample <= 1.0; N14_sample += step) {
			//System.out.println(C12_sample + "\t" + 0.0);
		}
		
		//System.out.println("Final Hit: " + simulated_N15_max_score + "\t" + max_score);
		double[] result = new double[2];
		result[0] = simulated_N15_max_score;
		result[1] = max_score;
		return result;
	}
	public static double[] simulate_c13(String formula, IsotopePattern c12_removed_clean_query_peaks_target, double ppm, int charge) {
		AIMPuritySampler sampler = new AIMPuritySampler();
		// find C13 peaks
		double big_step = 0.05;
		boolean doOnce = true;
		double smallest = 0;
		double largest = 1;
		double prevIndex = 0;
		double step = 0.01;
		double simulated_C13_max_score = -1;
		IsotopePatternSimilarity similarity = new IsotopePatternSimilarity();
		for (double C12_sample = 0; C12_sample <= 1; C12_sample += big_step) {
			IsotopePattern theoretical_isotope_c13 = sampler.calculate_pattern(formula, C12_sample, N14, H1, charge);
			IsotopePattern test_target = minimizeArea(c12_removed_clean_query_peaks_target, theoretical_isotope_c13, ppm);
			IsotopePattern test_target_norm = IsotopePatternManipulator.normalize(test_target);
			if (test_target_norm.getNumberOfIsotopes() > 0) {
				double score = similarity.compare(test_target_norm, theoretical_isotope_c13);
				if (score > 0) {
					if (doOnce) {
						smallest = prevIndex;
						doOnce = false;
					}
					largest = C12_sample;
				}
			}
			prevIndex = C12_sample;
		}
		largest = largest + big_step;
		if (largest <= 1 - big_step) {
			largest += big_step;
			if (largest >= 1) {
				largest = 1.0;
			}
		}
		double max_score = 0;
		for (double C12_sample = 0; C12_sample < smallest; C12_sample += step) {
			//System.out.println(C12_sample + "\t" + 0.0);
		}
		for (double C12_sample = smallest; C12_sample <= largest; C12_sample += step) {			
			IsotopePattern theoretical_isotope_c13 = sampler.calculate_pattern(formula, C12_sample, N14, H1, charge);
			IsotopePattern test_target = minimizeArea(c12_removed_clean_query_peaks_target, theoretical_isotope_c13, ppm);
			IsotopePattern test_target_norm = IsotopePatternManipulator.normalize(test_target);
			if (test_target_norm.getNumberOfIsotopes() > 0) {
				double score = similarity.compare(test_target_norm, theoretical_isotope_c13);
				if (max_score < score) {
					max_score = score;
					simulated_C13_max_score = C12_sample;
				}
				
				//System.out.println(C12_sample + "\t" + score);
			}
		}
		for (double C12_sample = largest + step; C12_sample <= 1.0; C12_sample += step) {
			//System.out.println(C12_sample + "\t" + 0.0);
		}
		
		//System.out.println("Final Hit: " + simulated_C13_max_score + "\t" + max_score);
		double[] result = new double[2];
		result[0] = simulated_C13_max_score;
		result[1] = max_score;
		return result;
	}
	public static IsotopePattern minimizeArea(IsotopePattern pattern, IsotopePattern theoretical_isotope, double ppm) {
		IsotopePattern result = new IsotopePattern();
		double minMass = Double.MAX_VALUE;
		double maxMass = Double.MIN_VALUE;								
		for (int i = 0; i < theoretical_isotope.getNumberOfIsotopes(); i++) {
			IsotopeContainer temp = theoretical_isotope.getIsotope(i);
			double mass = temp.getMass();
			double intensity = temp.getIntensity();
			if (intensity > 0.01) {
				double diff = (ppm * mass / 1E6);
				if (minMass > mass - diff) {
					minMass = mass - diff;
				}
				if (maxMass < mass + diff) {
					maxMass = mass + diff;
				}
			}
			//System.out.println(mass + "\t" + intensity);
		}
		result = minimizeArea(pattern, minMass, maxMass);
		for (int i = 0; i < result.getNumberOfIsotopes(); i++) {
			IsotopeContainer temp = result.getIsotope(i);
			double mass = temp.getMass();
			double intensity = temp.getIntensity();
			//System.out.println("Cleaned: " + mass + "\t" + intensity);
		}
		return result;
	}
	
	public static IsotopePattern removeC12Peak(IsotopePattern pattern, IsotopePattern C12_Peaks) {
		IsotopePattern result = new IsotopePattern();
		HashMap map = new HashMap();
		for (int j = 0; j < C12_Peaks.getNumberOfIsotopes(); j++) {
			IsotopeContainer container_C12 = C12_Peaks.getIsotope(j);
			double mass_c12 = container_C12.getMass();
			double intensity_c12 = container_C12.getIntensity();
			map.put(mass_c12, mass_c12);
		}
		for (int i = 0; i < pattern.getNumberOfIsotopes(); i++) {
			IsotopeContainer container = pattern.getIsotope(i);
			double mass = container.getMass();
			double intensity = container.getIntensity();
			if (!map.containsKey(mass)) {
				result.addIsotope(new IsotopeContainer(mass, intensity));
			}
			
		}
		return result;
	}
	public static IsotopePattern minimizeArea(IsotopePattern pattern, double minMass, double maxMass) {
		IsotopePattern result = new IsotopePattern();
		for (int i = 0; i < pattern.getNumberOfIsotopes(); i++) {
			IsotopeContainer container = pattern.getIsotope(i);
			double mass = container.getMass();
			double intensity = container.getIntensity();
			if (mass >= minMass && mass <= maxMass) {
				result.addIsotope(new IsotopeContainer(mass, intensity));
			}
		}
		return result;
	}
	public static IsotopePattern cleanPeak(IsotopePattern pattern, double distanceCutoff) {

		IsotopePattern result = new IsotopePattern();
		LinkedList list = new LinkedList();
		double prevMass = 0;
		double prevIntensity = 0;
		double highestPeakMass = 0;
		double highestPeakIntensity = 0;
		for (int i = 0; i < pattern.getNumberOfIsotopes(); i++) {			
			//System.out.println(peaks[0][i] + "\t" + peaks[1][i]);
			IsotopeContainer container = pattern.getIsotope(i);
			double mass = container.getMass();
			double intensity = container.getIntensity();
			if (mass - prevMass < distanceCutoff && intensity > prevIntensity) {
				if (list.size() > 0) {
					list.removeLast();
				}
				highestPeakMass = mass;
				highestPeakIntensity = intensity;
				list.add(mass);				
			} else if (mass - highestPeakMass < distanceCutoff && highestPeakIntensity > intensity){
				// skipping these
			} else {
				//if (intensity > 0) {
				list.add(mass);
				//}
			}
			
			prevMass = mass;
			prevIntensity = intensity;
			//								
		}
		//System.out.println("Reduced size: " + list.size());
		//System.out.println("Original size: " + pattern.getNumberOfIsotopes());
		for (int i = 0; i < pattern.getNumberOfIsotopes(); i++) {;			
			//System.out.println(peaks[0][i] + "\t" + peaks[1][i]);
			IsotopeContainer container = pattern.getIsotope(i);
			double mass = container.getMass();
			double intensity = container.getIntensity();
			if (list.contains(mass)) {		
				if (intensity > 0) {
					result.addIsotope(new IsotopeContainer(mass, intensity));
					//System.out.println("Cleaned: " + mass + "\t" + intensity);
				}
			}
		}
		return result;
	}
	
	public static IsotopePattern getPeakInfo(float[][] peaks ) {
		
		try {

			IsotopePattern result = new IsotopePattern();
			
			for (int i = 0; i < peaks[0].length; i++) {;			
				//System.out.println(peaks[0][i] + "\t" + peaks[1][i]);
				double mass = new Double(peaks[0][i]);
				double intensity = new Double(peaks[1][i]);
				result.addIsotope(new IsotopeContainer(mass, intensity));										
			}
			
			return result;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}*/
}

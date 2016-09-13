package PurityEstimation;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.LinkedList;

import org.openscience.cdk.formula.IsotopeContainer;
import org.openscience.cdk.formula.IsotopePattern;
import org.openscience.cdk.formula.IsotopePatternManipulator;

import CDKFunction.IsotopePatternSimilarity;
import METABOLOMICS.AIM.AIMPuritySampler;

/**
 * A toolset for matching peaks to theoretical peaks based on a chemical formula.
 * @author tshaw
 *
 */
public class MatchIsotopeToolbox {
	
	private static double C12 = 0.9893;
	private static double N14 = 0.99632;	
	private static double H1 = 0.999885;
	
	public static double[] simulate_n15(String formula, IsotopePattern c12_removed_clean_query_peaks_target, double ppm, int charge, double min, double max) {
		AIMPuritySampler sampler = new AIMPuritySampler();
		// find C13 peaks
		double big_step = 0.1;
		boolean doOnce = true;
		double smallest = 0;
		double largest = 1;
		double prevIndex = 0;
		double step = 0.05;
		double simulated_N15_max_score = -1;
		IsotopePatternSimilarity similarity = new IsotopePatternSimilarity();
		for (double N14_sample = min; N14_sample <= max; N14_sample += big_step) {
			IsotopePattern theoretical_isotope_n14 = sampler.calculate_pattern(formula, C12, N14_sample, H1, charge);
			IsotopePattern test_target = minimizeArea(c12_removed_clean_query_peaks_target, theoretical_isotope_n14, ppm);
			if (test_target.getNumberOfIsotopes() > 0) {
				try {
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
				} catch (Exception e) {
					//e.printStackTrace();
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
		if (min > smallest) {
			smallest = min;
		}
		if (max < largest) {
			largest = max;
		}
		for (double N14_sample = 0; N14_sample < smallest; N14_sample += step) {
			//System.out.println(C12_sample + "\t" + 0.0);
		}
		for (double N14_sample = smallest; N14_sample <= largest; N14_sample += step) {			
			IsotopePattern theoretical_isotope_c13 = sampler.calculate_pattern(formula, C12, N14_sample, H1, charge);
			IsotopePattern test_target = minimizeArea(c12_removed_clean_query_peaks_target, theoretical_isotope_c13, ppm);
			if (test_target.getNumberOfIsotopes() > 0) {
				try {
					IsotopePattern test_target_norm = IsotopePatternManipulator.normalize(test_target);
					if (test_target_norm.getNumberOfIsotopes() > 0) {
						double score = similarity.compare(test_target_norm, theoretical_isotope_c13);
						if (max_score < score) {
							max_score = score;
							simulated_N15_max_score = N14_sample;
						}
						
						//System.out.println(C12_sample + "\t" + score);
					}
				} catch (Exception e) {
					//e.printStackTrace();
				}
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
	public static double[] simulate_c13(String formula, IsotopePattern c12_removed_clean_query_peaks_target, double ppm, int charge, double min, double max) {
		AIMPuritySampler sampler = new AIMPuritySampler();
		// find C13 peaks
		double big_step = 0.1;
		boolean doOnce = true;
		double smallest = 0;
		double largest = 1;
		double prevIndex = 0;
		double step = 0.05;
		double simulated_C13_max_score = -1;
		IsotopePatternSimilarity similarity = new IsotopePatternSimilarity();
		for (double C12_sample = min; C12_sample <= max; C12_sample += big_step) {
			IsotopePattern theoretical_isotope_c13 = sampler.calculate_pattern(formula, C12_sample, N14, H1, charge);
			for (int i = 0; i < theoretical_isotope_c13.getNumberOfIsotopes(); i++) {
				IsotopeContainer temp = theoretical_isotope_c13.getIsotope(i);
				double mass = temp.getMass();
				double intensity = temp.getIntensity();
				if (intensity > 0.01) {
					//System.out.println("Theoretical: " + C12_sample + "\t" + mass + "\t" + intensity);
				}
			}
			IsotopePattern test_target = minimizeArea(c12_removed_clean_query_peaks_target, theoretical_isotope_c13, ppm);
			if (test_target.getNumberOfIsotopes() > 0) {
				try {
					IsotopePattern test_target_norm = IsotopePatternManipulator.normalize(test_target);
					for (int i = 0; i < test_target_norm.getNumberOfIsotopes(); i++) {
						IsotopeContainer temp = test_target_norm.getIsotope(i);
						double mass = temp.getMass();
						double intensity = temp.getIntensity();
						if (intensity > 0.01) {
							//System.out.println("C12_sample: " + C12_sample + "\t" + mass + "\t" + intensity);
						}
					}
					if (test_target_norm.getNumberOfIsotopes() > 0) {
						double score = similarity.compare(test_target_norm, theoretical_isotope_c13);
						if (score > 0) {
							if (doOnce) {
								smallest = prevIndex;
								doOnce = false;
							}
							largest = C12_sample;
							//System.out.println("Score C12_sample: " + C12_sample + "\t" + score);
						}
						
					}
				} catch (Exception e) {
					//e.printStackTrace();
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
		//for (double C12_sample = 0; C12_sample < smallest; C12_sample += step) {
		//	System.out.println(C12_sample + "\t" + 0.0);
		//}
		
		if (min > smallest) {
			smallest = min;
		}
		if (max < largest) {
			largest = max;
		}
		
		for (double C12_sample = smallest; C12_sample <= largest; C12_sample += step) {			
			IsotopePattern theoretical_isotope_c13 = sampler.calculate_pattern(formula, C12_sample, N14, H1, charge);
			IsotopePattern test_target = minimizeArea(c12_removed_clean_query_peaks_target, theoretical_isotope_c13, ppm);
			if (test_target.getNumberOfIsotopes() > 0) {
				try {
					IsotopePattern test_target_norm = IsotopePatternManipulator.normalize(test_target);
					if (test_target_norm.getNumberOfIsotopes() > 0) {
						double score = similarity.compare(test_target_norm, theoretical_isotope_c13);
						if (max_score < score) {
							max_score = score;
							simulated_C13_max_score = C12_sample;
						}
						
						//System.out.println(C12_sample + "\t" + score);
					}
				} catch (Exception e) {
					//e.printStackTrace();
				}
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
		//System.out.println("Min and Max Mass: " + minMass + "\t" + maxMass);
		IsotopePattern final_result = new IsotopePattern();
		result = minimizeAreaMinMax(pattern, minMass, maxMass);
		//System.out.println("After Minimizing number of isotopes: " + result.getNumberOfIsotopes());
		for (int i = 0; i < result.getNumberOfIsotopes(); i++) {
			IsotopeContainer temp = result.getIsotope(i);
			double mass = temp.getMass();
			double intensity = temp.getIntensity();
			//System.out.println("MinimizeArea: " + mass + "\t" + intensity);
			boolean find = false;
			for (int j = 0; j < theoretical_isotope.getNumberOfIsotopes(); j++) {
				IsotopeContainer theoretical_temp = theoretical_isotope.getIsotope(j);
				double theoretical_mass = theoretical_temp.getMass();
				double theoretical_intensity = theoretical_temp.getIntensity();
				double diff = (ppm * mass / 1E6);
				if (theoretical_mass >= mass - diff && theoretical_mass <= mass + diff && theoretical_intensity > 0.01) {
					find = true;
				}
			}
			if (find) {
				final_result.addIsotope(new IsotopeContainer(mass, intensity));
			}
			//System.out.println("Cleaned: " + mass + "\t" + intensity);
		}
		return final_result;
	}
	
	/**
	 * 
	 * @param pattern have more peaks
	 * @param C12_Peaks the target peaks to remove
	 * @return
	 */
	public static IsotopePattern removeC12Peak(IsotopePattern pattern, IsotopePattern C12_Peaks) {
		IsotopePattern result = new IsotopePattern();
		double max_mass = -1;
		double max_intensity = 0;
		double min_mass1 = -1;
		double min_mass2 = -1;
		double min_mass3 = -1;
		HashMap map = new HashMap();
		for (int j = 0; j < C12_Peaks.getNumberOfIsotopes(); j++) {
			IsotopeContainer container_C12 = C12_Peaks.getIsotope(j);
			double mass_c12 = container_C12.getMass();
			double intensity_c12 = container_C12.getIntensity();
			if (j == 0) {
				min_mass1 = mass_c12;
			} else if (j == 1) {
				min_mass2 = mass_c12;
			} 
			if (max_intensity < intensity_c12) {
				max_intensity = intensity_c12;
				max_mass = mass_c12;
			}
			//map.put(mass_c12, mass_c12);
		}
		map.put(max_mass, max_mass);
		for (int i = 0; i < pattern.getNumberOfIsotopes(); i++) {
			IsotopeContainer container = pattern.getIsotope(i);
			double mass = container.getMass();
			double intensity = container.getIntensity();
			if (!map.containsKey(mass) && mass != min_mass1 && mass != min_mass2) {
				result.addIsotope(new IsotopeContainer(mass, intensity));
			}
			
		}
		return result;
	}
	public static IsotopePattern minimizeAreaMinMax(IsotopePattern pattern, double minMass, double maxMass) {
		IsotopePattern result = new IsotopePattern();
		for (int i = 0; i < pattern.getNumberOfIsotopes(); i++) {
			
			IsotopeContainer container = pattern.getIsotope(i);
			double mass = container.getMass();
			double intensity = container.getIntensity();
			//
			if (mass >= minMass && mass <= maxMass) {
				//System.out.println("Check minimizeArea: " + mass + "\t" + intensity + "\tmin:" + minMass + "\tmax:" + maxMass);
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
	public static IsotopePattern readDTAFile(String inputFile) {
		try {
			IsotopePattern result = new IsotopePattern();
			
			FileInputStream fstream = new FileInputStream(inputFile);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));			
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				if (split.length >= 2) {					
					double mass = new Double(split[0]);
					double intensity = new Double(split[1]);
					result.addIsotope(new IsotopeContainer(mass, intensity));
				}
			}
			in.close();
			return result;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
	/**
	 * based on the input file, obtain the estimated mass
	 * @param fileName
	 * @return
	 */
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
	}
}

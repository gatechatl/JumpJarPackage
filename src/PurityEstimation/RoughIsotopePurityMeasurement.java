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
import MISC.ToolBox;

/**
 * Generate a rough isotope purity sampler
 * @author tshaw
 *
 */
public class RoughIsotopePurityMeasurement {
		
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
		return "[inputmzXMLFile] [scanNumber] [formula] [C12/C13/BOTH] [charge] [ppm] [minPurity] [maxPurity] [minMassCutoff] [maxMassCutoff]";
	}
	public static void main(String[] args) {
		System.out.println("mzXMLParser");
	}
	public static void execute(String[] args) {
		
		try {
			
			String inputFile = args[0];			
			int scanNumber = new Integer(args[1]);
			String formula = args[2]; // neutral form
			//double monoIsotopicMass = new Double(args[3]);
			String type = args[3];
			int charge = new Integer(args[4]);
			double ppm = new Double(args[5]);
			double minPurity = new Double(args[6]);
			double maxPurity = new Double(args[7]);
			double minMassCutoff = new Double(args[8]);
			double maxMassCutoff = new Double(args[9]);
			MSXMLParser parser = new MSXMLParser(inputFile);
			Scan scan = parser.rap(scanNumber);
			float[][] peaks = scan.getMassIntensityList();
			MZXMLFileInfo info = parser.getHeaderInfo();
			DataProcessingInfo dataProcessingInfo = info.getDataProcessing();
			IsotopePattern query_peaks = MatchIsotopeToolbox.getPeakInfo(peaks);
			long startTime = System.currentTimeMillis();
			/*for (int i = 0; i < query_peaks.getNumberOfIsotopes(); i++) {
				IsotopeContainer container = query_peaks.getIsotope(i);
				double mass = container.getMass();
				double intensity = container.getIntensity();
				System.out.println("Raw: " + mass + "\t" + intensity);
			}*/
			
			//IsotopePattern clean_query_peaks = MatchIsotopeToolbox.cleanPeak(query_peaks, 0.02); // clean the peaks
			IsotopePattern clean_query_peaks = MatchIsotopeToolbox.minimizeAreaMinMax(query_peaks, minMassCutoff, maxMassCutoff);
			//IsotopePattern clean_query_peaks = query_peaks;
			
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
			/*for (int i = 0; i < clean_query_peaks.getNumberOfIsotopes(); i++) {
				IsotopeContainer container = clean_query_peaks.getIsotope(i);
				double mass = container.getMass();
				double intensity = container.getIntensity();
				//System.out.println("cleanedPeaks\t" + mass + "\t" + intensity);
			}*/
			//System.out.println("### End reading peaks ###");
			double minMass = Double.MAX_VALUE;
			double maxMass = Double.MIN_VALUE;
			AIMPuritySampler sampler = new AIMPuritySampler();
						
			formula = addCharge2Formula(formula, charge);
			IsotopePattern theoretical_isotope = sampler.calculate_pattern(formula, C12, N14, H1, Math.abs(charge));
			
			//System.out.println("### End Generating Theoretical Isotope ###");
			
			IsotopePattern clean_query_peaks_target = MatchIsotopeToolbox.minimizeArea(clean_query_peaks, theoretical_isotope, ppm);
			IsotopePattern clean_query_peaks_target_norm = IsotopePatternManipulator.normalize(clean_query_peaks_target);
			//System.out.println("clean_query_peaks_target.getNumberOfIsotopes(): " + clean_query_peaks_target.getNumberOfIsotopes());
			//System.out.println("clean_query_peaks.getNumberOfIsotopes(): " + clean_query_peaks.getNumberOfIsotopes());
			IsotopePattern c12_removed_clean_query_peaks_target = MatchIsotopeToolbox.removeC12Peak(clean_query_peaks, clean_query_peaks_target);
			
			//IsotopeContainer container = theoretical_isotope.getIsotope(0);
			
			
			IsotopePatternSimilarity similarity = new IsotopePatternSimilarity();
			double similarity_score = -1;
			if (clean_query_peaks_target_norm.getNumberOfIsotopes() > 0) {
				similarity_score = similarity.compare(theoretical_isotope, clean_query_peaks_target_norm);
			}
			//System.out.println("Similarity_score: " + similarity_score);
			
			// find C13 peaks
			//System.out.println("### Searching C12 peak ###");
			//System.out.println("clean_query_peaks_target.getNumberOfIsotopes(): " + clean_query_peaks_target.getNumberOfIsotopes());
			if (minPurity < 0) {
				minPurity = 0;
			}
			if (maxPurity > 1) {
				maxPurity = 1;
			}
			double[] c12_simulation_result = MatchIsotopeToolbox.simulate_c13(formula, clean_query_peaks_target, ppm, charge, minPurity, maxPurity);			
			
			//System.out.println("### Searching C13 peak ###");
			//System.out.println("c12_removed_clean_query_peaks_target.getNumberOfIsotopes(): " + c12_removed_clean_query_peaks_target.getNumberOfIsotopes());
			double[] c13_simulation_result = MatchIsotopeToolbox.simulate_c13(formula, c12_removed_clean_query_peaks_target, ppm, charge, minPurity, maxPurity);
			
			//double[] n15_simulation_result = simulate_n15(formula, c12_removed_clean_query_peaks_target, ppm, charge);
			
			//System.out.println(formula + "\t" + similarity_score + "\t" + c12_simulation_result[0] + "\t" + c12_simulation_result[1] + "\t" + c13_simulation_result[0] + "\t" + c13_simulation_result[1]);
			if (type.equals("C12")) {
				System.out.println(formula + "\t" + charge + "\t" + type + "\tC12Purity:" + c12_simulation_result[0] + "\tMatchScore:" + c12_simulation_result[1]);
			} else if (type.equals("C13")) {
				System.out.println(formula + "\t" + charge + "\t" + type + "\tC12Purity:" + c13_simulation_result[0] + "\tMatchScore:" + c13_simulation_result[1]);
			} else if (type.equals("BOTH")){
				System.out.println(formula + "\t" + charge + "\t" + type + "\tC12Purity:" + c12_simulation_result[0] + "\tMatchScore:" + c12_simulation_result[1] + "\tC12Purity%:" + c13_simulation_result[0] + "\tMatchScore:" + c13_simulation_result[1]);
			} else {
				System.out.println(formula + "\t" + charge + "\t" + type + "\t" + similarity_score + "\t" + c12_simulation_result[0] + "\t" + c12_simulation_result[1] + "\t" + c13_simulation_result[0] + "\t" + c13_simulation_result[1]);
			}
			long endTime = System.currentTimeMillis();
			//System.out.println(endTime - startTime);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static String addCharge2Formula(String formula, int charge) {
		String newformula = formula;
		if (charge > 0) {
			for (int i = 0; i < charge; i++) {
				newformula = ToolBox.HillSystemOrder_ADD_H(newformula);
			}
		} else if (charge < 0) {
			for (int i = charge; i < 0; i++) {
				newformula = ToolBox.HillSystemOrder_Remove_H(newformula);
			}
		}
		return newformula;
	}
}


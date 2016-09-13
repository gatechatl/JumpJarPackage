package Specialized;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.SortedSet;
import java.util.TreeSet;

import org.openscience.cdk.formula.IsotopeContainer;
import org.openscience.cdk.formula.IsotopePattern;

import ISOTOPEDISTRIBUTION.IsotopeCalculator;
import ISOTOPEDISTRIBUTION.MS1.CalculateDegradedPeptideIsotopePattern;

public class JUMPfIsotopeDistribution {
	public static double C12 = 0.9893;
	public static double N14 = 0.99632;
	public static double H1 = 0.999885;

	public static void execute(String[] args) {
		
		try {
			String outputFile = args[1];
        	FileWriter fwriter = new FileWriter(outputFile);
            BufferedWriter out = new BufferedWriter(fwriter);
            
			String fileName = args[0];
			FileInputStream fstream = new FileInputStream(fileName);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				if (split.length > 9 && !split[0].equals("Peptides")) {
						
					String peptide = split[0].split("\\.")[1];
					String chemFormula = CalculateDegradedPeptideIsotopePattern.peptide2formula(peptide);
					IsotopePattern pattern = IsotopeCalculator.calculate_pattern(chemFormula, 50, C12, N14, H1, 1);
					
					HashMap isotope = new HashMap();
					for (int i = 0; i < pattern.getNumberOfIsotopes(); i++) {
						IsotopeContainer container = pattern.getIsotope(i);
						isotope.put(container.getMass(), container.getIntensity());													
					}
					
					double mass = CalculateDegradedPeptideIsotopePattern.getMassOfHighestPeak(pattern);
					String isotopes = "";
					boolean first = true;
					SortedSet<Double> keys = new TreeSet<Double>(isotope.keySet());
					int monoisotope = -1;
					int count = 0;
					for (double key : keys) {
						
						double value = (Double)isotope.get(key);
						if (mass == key) {
							monoisotope = count;
						}
						System.out.println(key + "\t" + value);
						if (first) {
							isotopes += key + ":" + value;
						} else {
							isotopes += "," + key + ":" + value;
						}
					    first = false;
					    count++;
					}
					out.write(str + "\t" + monoisotope + "\t" + isotopes + "\n");
					out.flush();
					System.out.println(peptide);
				}
			}
			in.close();
			out.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

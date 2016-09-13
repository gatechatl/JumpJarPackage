package Specialized;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;

public class CalculateMinMaxIsotopes {

	public static void main(String[] args) {
		
		try {
			
			double min_0 = Integer.MAX_VALUE;
			double max_0 = Integer.MIN_VALUE;
			double min_1 = Integer.MAX_VALUE;
			double max_1 = Integer.MIN_VALUE;
			double min_2 = Integer.MAX_VALUE;
			double max_2 = Integer.MIN_VALUE;
			double min_3 = Integer.MAX_VALUE;
			double max_3 = Integer.MIN_VALUE;
			double min_4 = Integer.MAX_VALUE;
			double max_4 = Integer.MIN_VALUE;
			double min_5 = Integer.MAX_VALUE;
			double max_5 = Integer.MIN_VALUE;
			double min_6 = Integer.MAX_VALUE;
			double max_6 = Integer.MIN_VALUE;
			
			String fileName = args[0];
			FileInputStream fstream = new FileInputStream(fileName);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				int mono = new Integer(split[2]);
				double mass = new Double(split[0]);
				if (mono == 0) {
					if (mass <= min_0) {
						min_0 = mass;
					}
					if (mass >= max_0) {
						max_0 = mass;
					}
				}
				if (mono == 1) {
					if (mass <= min_1) {
						min_1 = mass;
					}
					if (mass >= max_1) {
						max_1 = mass;
					}
				}
				if (mono == 2) {
					if (mass <= min_2) {
						min_2 = mass;
					}
					if (mass >= max_2) {
						max_2 = mass;
					}
				}
				if (mono == 3) {
					if (mass <= min_3) {
						min_3 = mass;
					}
					if (mass >= max_3) {
						max_3 = mass;
					}
				}
				if (mono == 4) {
					if (mass <= min_4) {
						min_4 = mass;
					}
					if (mass >= max_4) {
						max_4 = mass;
					}
				}
				if (mono == 5) {
					if (mass <= min_5) {
						min_5 = mass;
					}
					if (mass >= max_5) {
						max_5 = mass;
					}
				}
				if (mono == 6) {
					if (mass <= min_6) {
						min_6 = mass;
					}
					if (mass >= max_6) {
						max_6 = mass;
					}
				}
			}
			in.close();
			System.out.println("Monopeak = 0 range min = " + min_0 + " max = " + max_0);
			System.out.println("Monopeak = 1 range min = " + min_1 + " max = " + max_1);
			System.out.println("Monopeak = 2 range min = " + min_2 + " max = " + max_2);
			System.out.println("Monopeak = 3 range min = " + min_3 + " max = " + max_3);
			System.out.println("Monopeak = 4 range min = " + min_4 + " max = " + max_4);
			System.out.println("Monopeak = 5 range min = " + min_5 + " max = " + max_5);
			System.out.println("Monopeak = 6 range min = " + min_6 + " max = " + max_6);
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

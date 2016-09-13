package Specialized.XushengDrewPaper;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;

public class MISSILEGenerateFragmentationListScript {

	public static void execute(String[] args) {
		try {
			String inputFile = args[0];
			FileInputStream fstream = new FileInputStream(inputFile);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				String id = split[0];
				String smile = split[2];
				smile = smile.replaceAll("\\(", "\\\\(");
				smile = smile.replaceAll("\\)", "\\\\)");
				System.out.println("perl JumpJar.pl -program=SmileFragmenter -Smile='" + smile + "' -mass=50 -depth=2 > fragment_result_" + id);
			}
			in.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

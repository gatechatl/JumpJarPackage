package MZXMLParser;

import org.systemsbiology.jrap.DataProcessingInfo;
import org.systemsbiology.jrap.MSInstrumentInfo;
import org.systemsbiology.jrap.MSXMLParser;
import org.systemsbiology.jrap.MZXMLFileInfo;
import org.systemsbiology.jrap.Scan;

/**
 * mzXML parser into an object
 * @author tshaw
 *
 */
public class MzXMLParser {

	public static String description() {
		return "Parse a particular scan and generate the peak list";
	}
	public static String type() {
		return "MZXML";
	}
	public static String parameter_info() {
		return "[inputmzXMLFile] [scanNumber]";
	}
	public static void main(String[] args) {
		System.out.println("mzXMLParser");
	}
	public static void execute(String[] args) {
		
		try {
			
			String inputFile = args[0];
			int scanNumber = new Integer(args[1]);
			MSXMLParser parser = new MSXMLParser(inputFile);
			Scan scan = parser.rap(scanNumber);
			float[][] peaks = scan.getMassIntensityList();
			MZXMLFileInfo info = parser.getHeaderInfo();
			DataProcessingInfo dataProcessingInfo = info.getDataProcessing();
			
			MSInstrumentInfo instrumentInfo = info.getInstrumentInfo();
			System.out.println("Instrument Info: " + instrumentInfo.getModel());
			System.out.println(peaks[0].length);
			System.out.println(peaks[1].length);
			for (int i = 0; i < peaks[0].length; i++) {;			
				System.out.println(peaks[0][i] + "\t" + peaks[1][i]);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

package MassDatabaseGeneration;

public class CreateIndexFormulaOriginalScript {

	public static void execute(String[] args) {
		
		String inputFile = args[0];
		String outputFolder = args[1];
		int parallel_num = new Integer(args[2]);
		int buffer = 1600 / parallel_num;
		for (int i = 0; i <= 1600; i = i + buffer) {
			int max = i + buffer;
			if (max > 1600) {
				max = 1600;
			}
			System.out.println("jumpjar -CreateIndexFormulaOriginal " + inputFile + " " + outputFolder + " " + i + " " + (i + buffer));
		}
	}
}

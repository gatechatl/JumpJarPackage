package METABOLOMICS.AIM;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.util.List;

import Interface.JumpInterface;

public class AIMFragmenter {
	
	public static void main(String[] args) {

		//String smile = "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)OC4C(C(C(C(O4)CO)O)O)O)O)O";
		String smile = "C(CC(=O)NC(CS)C(=O)NCC(=O)O)C(C(=O)[O-])[NH3+]";
		double mass = 18;
		try {
			MetFragFragmenter(smile, mass, 1, false, true);
		
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void MetFragFragmenter(String smile, double mass, int depth, boolean aromatic_ring_flag, boolean molecularFormulaRedundancyCheck) {
		FragmentSingleCompound fragmenter = new FragmentSingleCompound();
		fragmenter.calculate(smile, mass, depth, aromatic_ring_flag, molecularFormulaRedundancyCheck);		
	}
	
	public static String parameter_info() {
		return "[smile string] [mass (cutoff)] [depth (number)] [aromatic_ring_flag yes/no] [neutralLossFile] [bondEnergyFile]";
	}
	public static void execute(String[] args) {
		String smile = args[0]; //"C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)OC4C(C(C(C(O4)CO)O)O)O)O)O";
		
		double mass = new Double(args[1]); //303.0499;
		
		int depth = new Integer(args[2]);
		
		boolean aromatic_ring_flag = false;
		String aromatic_ring_str = args[3];
		if (aromatic_ring_str.equals("yes")) {
			aromatic_ring_flag = true;
		}
		String molecularFormulaRedundancyCheck_str = args[4];
		boolean molecularFormulaRedundancyCheck = true;
		if (molecularFormulaRedundancyCheck_str.equals("no")) {
			molecularFormulaRedundancyCheck = false;
		}
		String neutralLossFile = args[5];	
		String bondEnergyFile = args[6];
		
		File neutral = new File(neutralLossFile);
		File new_netrual = new File("neutralLoss.csv");
		if (!new_netrual.exists() && neutral.exists()) {
			try {
				//System.out.println("Copied file to neutralLoss.csv");
				copyFile(neutral, new_netrual);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				System.out.println("Failed to copy");
			}
		} else {
			
			//System.out.println("Can't find Neutral Loss File: " + neutralLossFile);
		}
		File bond = new File(bondEnergyFile);
		File new_bond = new File("bondenergies.txt");
		if (!new_bond.exists() && bond.exists()) {
			try {
				copyFile(bond, new_bond);
				//System.out.println("Copied file to bondenergies.txt");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				//System.out.println("Failed to copy");
			}
		} else {
			//System.out.println("Can't find Bond energy File: " + bondEnergyFile);
		}
		
		MetFragFragmenter(smile, mass, depth, aromatic_ring_flag, molecularFormulaRedundancyCheck);
	}
	
	public static void copyFile(File sourceFile, File destFile) throws IOException {
	    if(!destFile.exists()) {
	        destFile.createNewFile();
	    }

	    FileChannel source = null;
	    FileChannel destination = null;

	    try {
	        source = new FileInputStream(sourceFile).getChannel();
	        destination = new FileOutputStream(destFile).getChannel();
	        destination.transferFrom(source, 0, source.size());
	    }
	    finally {
	        if(source != null) {
	            source.close();
	        }
	        if(destination != null) {
	            destination.close();
	        }
	    }
	}
}

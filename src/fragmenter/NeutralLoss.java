package fragmenter;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.interfaces.IMolecularFormula;

public class NeutralLoss {
	
	private IMolecularFormula elementalComposition;
	private int mode;
	private IMolecularFormula topoFragment;
	private int hydrogenDifference;
	private int distance;
	private String atomToStart;
	private int hydrogenOnStartAtom;
	private List<String> mainAtoms;
		
	/**
	 * Instantiates a new neutral loss.
	 * Mode
	 * 
	 * @param elementalComposition the elemental composition
	 * @param mode the mode
	 * @param topoFragment the topo fragment
	 * @param hydrogenDifference the hydrogen difference
	 * @param distance the distance
	 * @param atomToStart the atom to start
	 * @param hydrogenOnStartAtom the hydrogen on start atom
	 */
	public NeutralLoss(IMolecularFormula elementalComposition, IMolecularFormula topoFragment, int mode, int hydrogenDifference, int distance, String atomToStart, int hydrogenOnStartAtom)
	{
		setElementalComposition(elementalComposition);
		setMode(mode); 
		setTopoFragment(topoFragment);
		setHydrogenDifference(hydrogenDifference);
		setDistance(distance);
		setAtomToStart(atomToStart);
		setHydrogenOnStartAtom(hydrogenOnStartAtom);
		mainAtomSet(topoFragment);
	}
	
	
	/**
	 * Find the main atoms in the neutral loss
	 * 
	 * @param topoFragment the topo fragment
	 */
	private void mainAtomSet(IMolecularFormula topoFragment)
	{
		mainAtoms = new ArrayList<String>();
		for (IIsotope isotope : topoFragment.isotopes()) {
			int count = this.topoFragment.getIsotopeCount(isotope);
    		for (int i = 0; i < count; i++)
    		{
    			if(!isotope.getSymbol().equals("H"))
    				mainAtoms.add(isotope.getSymbol());
    		}
		}
	}

	/**
	 * Sets the elemental composition.
	 * 
	 * @param elementalComposition the new elemental composition
	 */
	public void setElementalComposition(IMolecularFormula elementalComposition) {
		this.elementalComposition = elementalComposition;
	}

	/**
	 * Gets the elemental composition.
	 * 
	 * @return the elemental composition
	 */
	public IMolecularFormula getElementalComposition() {
		return elementalComposition;
	}

	/**
	 * Sets the mode.
	 * 
	 * @param mode the new mode
	 */
	public void setMode(int mode) {
		this.mode = mode;
	}

	/**
	 * Gets the mode.
	 * 
	 * @return the mode
	 */
	public int getMode() {
		return mode;
	}

	public void setHydrogenDifference(int hydrogenDifference) {
		this.hydrogenDifference = hydrogenDifference;
	}

	public int getHydrogenDifference() {
		return hydrogenDifference;
	}

	public void setTopoFragment(IMolecularFormula topoFragment) {
		this.topoFragment = topoFragment;
	}

	public IMolecularFormula getTopoFragment() {
		return topoFragment;
	}

	public void setDistance(int distance) {
		this.distance = distance;
	}

	public int getDistance() {
		return distance;
	}

	public void setAtomToStart(String atomToStart) {
		this.atomToStart = atomToStart;
	}

	public String getAtomToStart() {
		return atomToStart;
	}

	public void setHydrogenOnStartAtom(int hydrogenOnStartAtom) {
		this.hydrogenOnStartAtom = hydrogenOnStartAtom;
	}

	public int getHydrogenOnStartAtom() {
		return hydrogenOnStartAtom;
	}

	public void setMainAtoms(List<String> mainAtoms) {
		this.mainAtoms = mainAtoms;
	}

	public List<String> getMainAtoms() {
		return mainAtoms;
	}

}

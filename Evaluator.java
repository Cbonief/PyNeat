package com.hydrozoa.hydroneat;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;

public abstract class Evaluator {
	
	private FitnessGenomeComparator fitComp = new FitnessGenomeComparator();
	
	private Counter nodeInnovation;
	private Counter connectionInnovation;
	
	private Random random = new Random();
	
	/* Constants for tuning */
	private float C1 = 1.0f;
	private float C2 = 1.0f;
	private float C3 = 0.4f;
	private float DT = 10.0f;
	private float MUTATION_RATE = 0.5f;
	private float ADD_CONNECTION_RATE = 0.1f;
	private float ADD_NODE_RATE = 0.1f;
	
	private int populationSize;
	
	private List<Genome> genomes;
	private List<Genome> nextGenGenomes;
	
	private List<Species> species;
	
	private Map<Genome, Species> mappedSpecies;
	private Map<Genome, Float> scoreMap;
	private float highestScore;
	private Genome fittestGenome;
	
	public Evaluator(int populationSize, Genome startingGenome, Counter nodeInnovation, Counter connectionInnovation) {
		this.populationSize = populationSize;
		this.nodeInnovation = nodeInnovation;
		this.connectionInnovation = connectionInnovation;
		genomes = new ArrayList<Genome>(populationSize);
		for (int i = 0; i < populationSize; i++) {
			genomes.add(new Genome(startingGenome));
		}
		nextGenGenomes = new ArrayList<Genome>(populationSize);
		mappedSpecies = new HashMap<Genome, Species>();
		scoreMap = new HashMap<Genome, Float>();
		species = new ArrayList<Species>();
	}
	
	/**
	 * Runs one generation
	 */
	public void evaluate() {
		// Reset everything for next generation
		for (Species s : species) {
			s.reset(random);
		}
		scoreMap.clear();
		mappedSpecies.clear();
		nextGenGenomes.clear();
		highestScore = Float.MIN_VALUE;
		fittestGenome = null;
		
		// Place genomes into species
		for (Genome g : genomes) {
			boolean foundSpecies = false;
			for (Species s : species) {
				if (Genome.compatibilityDistance(g, s.mascot, C1, C2, C3) < DT) { // compatibility distance is less than DT, so genome belongs to this species
					s.members.add(g);
					mappedSpecies.put(g, s);
					foundSpecies = true;
					break;
				}
			}
			if (!foundSpecies) { // if there is no appropiate species for genome, make a new one
				Species newSpecies = new Species(g);
				species.add(newSpecies);
				mappedSpecies.put(g, newSpecies);
			}
		}
		
		// Remove unused species
		Iterator<Species> iter = species.iterator();
		while(iter.hasNext()) {
			Species s = iter.next();
			if (s.members.isEmpty()) {
				iter.remove();
			}
		}
		
		// Evaluate genomes and assign score
		for (Genome g : genomes) {
			Species s = mappedSpecies.get(g);		// Get species of the genome
			
			float score = evaluateGenome(g);
			float adjustedScore = score / mappedSpecies.get(g).members.size();
			
			s.addAdjustedFitness(adjustedScore);	
			s.fitnessPop.add(new FitnessGenome(g, adjustedScore));
			scoreMap.put(g, adjustedScore);
			if (score > highestScore) {
				highestScore = score;
				fittestGenome = g;
			}
		}
		
		// put best genomes from each species into next generation
		for (Species s : species) {
			Collections.sort(s.fitnessPop, fitComp);
			Collections.reverse(s.fitnessPop);
			FitnessGenome fittestInSpecies = s.fitnessPop.get(0);
			nextGenGenomes.add(fittestInSpecies.genome);
		}
		
		// Breed the rest of the genomes
		while (nextGenGenomes.size() < populationSize) { // replace removed genomes by randomly breeding
			Species s = getRandomSpeciesBiasedAjdustedFitness(random);
			
			Genome p1 = getRandomGenomeBiasedAdjustedFitness(s, random);
			Genome p2 = getRandomGenomeBiasedAdjustedFitness(s, random);
			
			Genome child;
			if (scoreMap.get(p1) >= scoreMap.get(p2)) {
				child = Genome.crossover(p1, p2, random);
			} else {
				child = Genome.crossover(p2, p1, random);
			}
			if (random.nextFloat() < MUTATION_RATE) {
				child.mutation(random);
			}
			if (random.nextFloat() < ADD_CONNECTION_RATE) {
				//System.out.println("Adding connection mutation...");
				child.addConnectionMutation(random, connectionInnovation, 10);
			}
			if (random.nextFloat() < ADD_NODE_RATE) {
				//System.out.println("Adding node mutation...");
				child.addNodeMutation(random, connectionInnovation, nodeInnovation);
			}
			nextGenGenomes.add(child);
		}
		
		genomes = nextGenGenomes;
		nextGenGenomes = new ArrayList<Genome>();
	}
	
	/**
	 * Selects a random species from the species list, where species with a higher total adjusted fitness have a higher chance of being selected
	 */
	private Species getRandomSpeciesBiasedAjdustedFitness(Random random) {
		double completeWeight = 0.0;	// sum of probablities of selecting each species - selection is more probable for species with higher fitness
		for (Species s : species) {
            completeWeight += s.totalAdjustedFitness;
		}
        double r = Math.random() * completeWeight;
        double countWeight = 0.0;
        for (Species s : species) {
            countWeight += s.totalAdjustedFitness;
            if (countWeight >= r) {
            	 return s;
            }
        }
        throw new RuntimeException("Couldn't find a species... Number is species in total is "+species.size()+", and the total adjusted fitness is "+completeWeight);
	}
	
	/**
	 * Selects a random genome from the species chosen, where genomes with a higher adjusted fitness have a higher chance of being selected
	 */
	private Genome getRandomGenomeBiasedAdjustedFitness(Species selectFrom, Random random) {
		double completeWeight = 0.0;	// sum of probablities of selecting each genome - selection is more probable for genomes with higher fitness
		for (FitnessGenome fg : selectFrom.fitnessPop) {
			completeWeight += fg.fitness;
		}
        double r = Math.random() * completeWeight;
        double countWeight = 0.0;
        for (FitnessGenome fg : selectFrom.fitnessPop) {
            countWeight += fg.fitness;
            if (countWeight >= r) {
            	 return fg.genome;
            }
        }
        throw new RuntimeException("Couldn't find a genome... Number is genomes in selæected species is "+selectFrom.fitnessPop.size()+", and the total adjusted fitness is "+completeWeight);
	}
	
	public int getSpeciesAmount() {
		return species.size();
	}
	
	public float getHighestFitness() {
		return highestScore;
	}
	
	public Genome getFittestGenome() {
		return fittestGenome;
	}
	
	protected abstract float evaluateGenome(Genome genome);
	
	public class FitnessGenome {
		
		float fitness;
		Genome genome;
		
		public FitnessGenome(Genome genome, float fitness) {
			this.genome = genome;
			this.fitness = fitness;
		}
	}
	
	public class Species {
		
		public Genome mascot;
		public List<Genome> members;
		public List<FitnessGenome> fitnessPop;
		public float totalAdjustedFitness = 0f;
		
		public Species(Genome mascot) {
			this.mascot = mascot;
			this.members = new LinkedList<Genome>(); 
			this.members.add(mascot);
			this.fitnessPop = new ArrayList<FitnessGenome>(); 
		}
		
		public void addAdjustedFitness(float adjustedFitness) {
			this.totalAdjustedFitness += adjustedFitness;
		}
		
		/*
		 *	 Selects new random mascot + clear members + set totaladjustedfitness to 0f
		 */
		public void reset(Random r) {
			int newMascotIndex = r.nextInt(members.size());
			this.mascot = members.get(newMascotIndex);
			members.clear();
			fitnessPop.clear();
			totalAdjustedFitness = 0f;
		}
	}
	
	public class FitnessGenomeComparator implements Comparator<FitnessGenome> {

		@Override
		public int compare(FitnessGenome one, FitnessGenome two) {
			if (one.fitness > two.fitness) {
				return 1;
			} else if (one.fitness < two.fitness) {
				return -1;
			}
			return 0;
		}
		
	}
}

package com.hydrozoa.hydroneat;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import com.hydrozoa.hydroneat.NodeGene.TYPE;

public class Genome {
	
	private static List<Integer> tmpList1 = new ArrayList<Integer>();
	private static List<Integer> tmpList2 = new ArrayList<Integer>();
	
	private final float PROBABILITY_PERTURBING = 0.9f; // rest is probability of assigning new weight
	
	private Map<Integer, ConnectionGene> connections;
	private Map<Integer, NodeGene> nodes;
	
	public Genome() {
		nodes = new HashMap<Integer, NodeGene>();
		connections = new HashMap<Integer, ConnectionGene>();
	}
	
	public Genome(Genome toBeCopied) {
		nodes = new HashMap<Integer, NodeGene>();
		connections = new HashMap<Integer, ConnectionGene>();
		
		for (Integer index : toBeCopied.getNodeGenes().keySet()) {
			nodes.put(index, new NodeGene(toBeCopied.getNodeGenes().get(index)));
		}
		
		for (Integer index : toBeCopied.getConnectionGenes().keySet()) {
			connections.put(index, new ConnectionGene(toBeCopied.getConnectionGenes().get(index)));
		}
	}
	
	public void addNodeGene(NodeGene gene) {
		nodes.put(gene.getId(), gene);
	}
	
	public void addConnectionGene(ConnectionGene gene) {
		connections.put(gene.getInnovation(), gene);
	}
	
	public Map<Integer, ConnectionGene> getConnectionGenes() {
		return connections;
	}
	
	public Map<Integer, NodeGene> getNodeGenes() {
		return nodes;
	}
	
	public void mutation(Random r) {
		for(ConnectionGene con : connections.values()) {
			if (r.nextFloat() < PROBABILITY_PERTURBING) { 			// uniformly perturbing weights
				con.setWeight(con.getWeight()*(r.nextFloat()*4f-2f));
			} else { 												// assigning new weight
				con.setWeight(r.nextFloat()*4f-2f);
			}
		}
	}
	
	public void addConnectionMutation(Random r, Counter innovation, int maxAttempts) {
		int tries = 0;
		boolean success = false;
		while (tries < maxAttempts && success == false) {
			tries++;

			// Cria um vetor que armazena os números de inovação dos nós.
			Integer[] nodeInnovationNumbers = new Integer[nodes.keySet().size()];
			nodes.keySet().toArray(nodeInnovationNumbers);
			Integer keyNode1 = nodeInnovationNumbers[r.nextInt(nodeInnovationNumbers.length)];
			Integer keyNode2 = nodeInnovationNumbers[r.nextInt(nodeInnovationNumbers.length)];
			
			NodeGene node1 = nodes.get(keyNode1);
			NodeGene node2 = nodes.get(keyNode2);
			float weight = r.nextFloat()*2f-1f;
			
			boolean reversed = false;
			if (node1.getType() == NodeGene.TYPE.HIDDEN && node2.getType() == NodeGene.TYPE.INPUT) {
				reversed = true;
			} else if (node1.getType() == NodeGene.TYPE.OUTPUT && node2.getType() == NodeGene.TYPE.HIDDEN) {
				reversed = true;
			} else if (node1.getType() == NodeGene.TYPE.OUTPUT && node2.getType() == NodeGene.TYPE.INPUT) {
				reversed = true;
			}
			
			boolean connectionImpossible = false;
			if (node1.getType() == NodeGene.TYPE.INPUT && node2.getType() == NodeGene.TYPE.INPUT) {
				connectionImpossible = true;
			} else if (node1.getType() == NodeGene.TYPE.OUTPUT && node2.getType() == NodeGene.TYPE.OUTPUT) {
				connectionImpossible = true;
			}
			
			boolean connectionExists = false;
			for (ConnectionGene con : connections.values()) {
				if (con.getInNode() == node1.getId() && con.getOutNode() == node2.getId()) { // existing connection
					connectionExists = true;
					break;
				} else if (con.getInNode() == node2.getId() && con.getOutNode() == node1.getId()) { // existing connection
					connectionExists = true;
					break;
				}
			}
			
			if (connectionExists || connectionImpossible) {
				continue;
			}
			
			ConnectionGene newCon = new ConnectionGene(reversed ? node2.getId() : node1.getId(), reversed ? node1.getId() : node2.getId(), weight, true, innovation.getInnovation());
			connections.put(newCon.getInnovation(), newCon);
			success = true;
		}
		if (success == false) {
			System.out.println("Tried, but could not add more connections");
		}
	}
	
	public void addNodeMutation(Random r, Counter connectionInnovation, Counter nodeInnovation) {
		ConnectionGene con = (ConnectionGene) connections.values().toArray()[r.nextInt(connections.size())];
		
		NodeGene inNode = nodes.get(con.getInNode());
		NodeGene outNode = nodes.get(con.getOutNode());
		
		con.disable();
		
		NodeGene newNode = new NodeGene(TYPE.HIDDEN, nodeInnovation.getInnovation());
		ConnectionGene inToNew = new ConnectionGene(inNode.getId(), newNode.getId(), 1f, true, connectionInnovation.getInnovation());
		ConnectionGene newToOut = new ConnectionGene(newNode.getId(), outNode.getId(), con.getWeight(), true, connectionInnovation.getInnovation());
		
		nodes.put(newNode.getId(), newNode);
		connections.put(inToNew.getInnovation(), inToNew);
		connections.put(newToOut.getInnovation(), newToOut);
	}
	
	/**
	 * @param parent1	More fit parent
	 * @param parent2	Less fit parent
	 */
	public static Genome crossover(Genome parent1, Genome parent2, Random r) {
		Genome child = new Genome();
		
		for (NodeGene parent1Node : parent1.getNodeGenes().values()) {
			child.addNodeGene(new NodeGene(parent1Node));
		}
		
		for (ConnectionGene parent1Node : parent1.getConnectionGenes().values()) {
			if (parent2.getConnectionGenes().containsKey(parent1Node.getInnovation())) { // matching gene
				ConnectionGene childConGene = r.nextBoolean() ? new ConnectionGene(parent1Node) : new ConnectionGene(parent2.getConnectionGenes().get(parent1Node.getInnovation()));
				child.addConnectionGene(childConGene);
			} else { // disjoint or excess gene
				ConnectionGene childConGene = new ConnectionGene(parent1Node);
				child.addConnectionGene(childConGene);
			}
		}
		
		return child;
	}
	
	public static float compatibilityDistance(Genome genome1, Genome genome2, float c1, float c2, float c3) {
		int excessGenes = countExcessGenes(genome1, genome2);
		int disjointGenes = countDisjointGenes(genome1, genome2);
		float avgWeightDiff = averageWeightDiff(genome1, genome2);
		
		return excessGenes * c1 + disjointGenes * c2 + avgWeightDiff * c3;
	}
	
	public static int countMatchingGenes(Genome genome1, Genome genome2) {
		int matchingGenes = 0;
		
		List<Integer> nodeKeys1 = asSortedList(genome1.getNodeGenes().keySet(), tmpList1);
		List<Integer> nodeKeys2 = asSortedList(genome2.getNodeGenes().keySet(), tmpList2);
		
		int highestInnovation1 = nodeKeys1.get(nodeKeys1.size()-1);
		int highestInnovation2 = nodeKeys2.get(nodeKeys2.size()-1);
		int indices = Math.max(highestInnovation1, highestInnovation2);
		
		for (int i = 0; i <= indices; i++) { 					// loop through genes -> i is innovation numbers
			NodeGene node1 = genome1.getNodeGenes().get(i);
			NodeGene node2 = genome2.getNodeGenes().get(i);
			if (node1 != null && node2 != null) { 
				// both genomes has the gene w/ this innovation number
				matchingGenes++;
			}
		}
		
		List<Integer> conKeys1 = asSortedList(genome1.getConnectionGenes().keySet(), tmpList1);
		List<Integer> conKeys2 = asSortedList(genome2.getConnectionGenes().keySet(), tmpList2);
		
		highestInnovation1 = conKeys1.get(conKeys1.size()-1);
		highestInnovation2 = conKeys2.get(conKeys2.size()-1);
		
		indices = Math.max(highestInnovation1, highestInnovation2);
		for (int i = 0; i <= indices; i++) { 					// loop through genes -> i is innovation numbers
			ConnectionGene connection1 = genome1.getConnectionGenes().get(i);
			ConnectionGene connection2 = genome2.getConnectionGenes().get(i);
			if (connection1 != null && connection2 != null) { 
				// both genomes has the gene w/ this innovation number
				matchingGenes++;
			}
		}
		
		return matchingGenes;
	}
	
	public static int countDisjointGenes(Genome genome1, Genome genome2) {
		int disjointGenes = 0;
		
		List<Integer> nodeKeys1 = asSortedList(genome1.getNodeGenes().keySet(), tmpList1);
		List<Integer> nodeKeys2 = asSortedList(genome2.getNodeGenes().keySet(), tmpList2);
		
		int highestInnovation1 = nodeKeys1.get(nodeKeys1.size()-1);
		int highestInnovation2 = nodeKeys2.get(nodeKeys2.size()-1);
		int indices = Math.max(highestInnovation1, highestInnovation2);
		
		for (int i = 0; i <= indices; i++) {
			NodeGene node1 = genome1.getNodeGenes().get(i);
			NodeGene node2 = genome2.getNodeGenes().get(i);
			if (node1 == null && highestInnovation1 > i && node2 != null) {
				// genome 1 lacks gene, genome 2 has gene, genome 1 has more genes w/ higher innovation numbers
				disjointGenes++;
			} else if (node2 == null && highestInnovation2 > i && node1 != null) {
				disjointGenes++;
			}
		}
		
		List<Integer> conKeys1 = asSortedList(genome1.getConnectionGenes().keySet(), tmpList1);
		List<Integer> conKeys2 = asSortedList(genome2.getConnectionGenes().keySet(), tmpList2);
		
		highestInnovation1 = conKeys1.get(conKeys1.size()-1);
		highestInnovation2 = conKeys2.get(conKeys2.size()-1);
		
		indices = Math.max(highestInnovation1, highestInnovation2);
		for (int i = 0; i <= indices; i++) {
			ConnectionGene connection1 = genome1.getConnectionGenes().get(i);
			ConnectionGene connection2 = genome2.getConnectionGenes().get(i);
			if (connection1 == null && highestInnovation1 > i && connection2 != null) {
				disjointGenes++;
			} else if (connection2 == null && highestInnovation2 > i && connection1 != null) {
				disjointGenes++;
			}
		}
		
		return disjointGenes;
	}
	
	public static int countExcessGenes(Genome genome1, Genome genome2) {
		int excessGenes = 0;
		
		List<Integer> nodeKeys1 = asSortedList(genome1.getNodeGenes().keySet(), tmpList1);
		List<Integer> nodeKeys2 = asSortedList(genome2.getNodeGenes().keySet(), tmpList2);
		
		int highestInnovation1 = nodeKeys1.get(nodeKeys1.size()-1);
		int highestInnovation2 = nodeKeys2.get(nodeKeys2.size()-1);
		int indices = Math.max(highestInnovation1, highestInnovation2);
		
		for (int i = 0; i <= indices; i++) {
			NodeGene node1 = genome1.getNodeGenes().get(i);
			NodeGene node2 = genome2.getNodeGenes().get(i);
			if (node1 == null && highestInnovation1 < i && node2 != null) {
				excessGenes++;
			} else if (node2 == null && highestInnovation2 < i && node1 != null) {
				excessGenes++;
			}
		}
		
		List<Integer> conKeys1 = asSortedList(genome1.getConnectionGenes().keySet(), tmpList1);
		List<Integer> conKeys2 = asSortedList(genome2.getConnectionGenes().keySet(), tmpList2);
		
		highestInnovation1 = conKeys1.get(conKeys1.size()-1);
		highestInnovation2 = conKeys2.get(conKeys2.size()-1);
		
		indices = Math.max(highestInnovation1, highestInnovation2);
		for (int i = 0; i <= indices; i++) {
			ConnectionGene connection1 = genome1.getConnectionGenes().get(i);
			ConnectionGene connection2 = genome2.getConnectionGenes().get(i);
			if (connection1 == null && highestInnovation1 < i && connection2 != null) {
				excessGenes++;
			} else if (connection2 == null && highestInnovation2 < i && connection1 != null) {
				excessGenes++;
			}
		}
		
		return excessGenes;
	}
	
	public static float averageWeightDiff(Genome genome1, Genome genome2) {
		int matchingGenes = 0;
		float weightDifference = 0;
		
		List<Integer> conKeys1 = asSortedList(genome1.getConnectionGenes().keySet(), tmpList1);
		List<Integer> conKeys2 = asSortedList(genome2.getConnectionGenes().keySet(), tmpList2);
		
		int highestInnovation1 = conKeys1.get(conKeys1.size()-1);
		int highestInnovation2 = conKeys2.get(conKeys2.size()-1);
		
		int indices = Math.max(highestInnovation1, highestInnovation2);
		for (int i = 0; i <= indices; i++) { 					// loop through genes -> i is innovation numbers
			ConnectionGene connection1 = genome1.getConnectionGenes().get(i);
			ConnectionGene connection2 = genome2.getConnectionGenes().get(i);
			if (connection1 != null && connection2 != null) { 
				// both genomes has the gene w/ this innovation number
				matchingGenes++;
				weightDifference += Math.abs(connection1.getWeight()-connection2.getWeight());
			}
		}
		
		return weightDifference/matchingGenes;
	}
	
	/**
	 * Note: Will sort in ascending order
	 */
	private static List<Integer> asSortedList(Collection<Integer> c, List<Integer> list) {
	  list.clear();
	  list.addAll(c);
	  java.util.Collections.sort(list);
	  return list;
	}
}

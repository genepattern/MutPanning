/************************************************************           
 * MutPanning - Main Class									*
 * 															*   
 * Author:		Felix Dietlein								*   
 *															*   
 * Copyright:	(C) 2018 									*   
 *															*   
 * License:		Public Domain								*   
 *															*   
 * Summary: MutPanning is detects cancer genes, based on	*
 * their mutation counts as well as the sequence context 	*
 * around mutations. This script coordinates the execution	*
 * of all substeps. All arguments are facultative except	*
 * of the root file in which the script should write its 	*
 * (intermediate) output files. In this folde the script	*
 * expects a "Hg19" folder and a "scr" folder. Please 		*
 * launch this script with -Xmx55G in the VM arguments to	*
 * allocate enough memory. Designed to run on a Linux		*
 * platform with at least 24 CPUs, 64 GB RAM and 2TB hard	*
 * disk space to store the intermediate output files.		*
 *************************************************************/



import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.lang.management.ManagementFactory;
import java.lang.management.OperatingSystemMXBean;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.Collections;

import com.sun.management.OperatingSystemMXBean.*;

public class MutPanning {
	
	/*	argument 0: root file, where all the other files can be found
	 * argument 1: maf file (standard value: root file/MutationsComplete.maf)
	 * argument 2: sample annotation file (standard value: root file/SamplesComplete.txt)
	 * argument3 : minimal no. samples per cluster (standard value 3)
	 * argument4: minimal no. mutations per cluster (standard value 1000)
	 *argument5: no cpus available to distribute clustering & CBASE workload (standard value 24)
	 arugment 6: min no. samples for CBASE (standard value 100)
	 argument 7: min no. mutations for CBASE (standard value 5000)
	argument 8: path to python script to compute distribution of nonsynonymous mutation counts (standard value: root file/src/ComputeDistribution.py)
	argument 9: path to Hg19 folder (standard value root file/Hg19/)
	
	 */
	
	static String[] index_header_maf={"Hugo_Symbol","Chromosome","Start_Position","End_Position","Strand","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2","Tumor_Sample_Barcode"};
	static String[] index_header_samples={"ID","Sample","Cohort"};
	
	//TODO: make the version variable so that not static columns are expected but rather by variable column headers
	//that we can write into the manual
	
	public static void main(String[] arguments){
		
		
		
		//fill arguments with standard values if not parsed. Only the root file cannot be replaced.
		String[] args=new String[11];
		for (int i=0;i<args.length;i++){
			args[i]="";
		}
		for (int i=0;i<arguments.length;i++){
			if(!arguments[i].equals("")&&!arguments[i].equals(" ")){
				args[i]=arguments[i];
			}
		}
		
		
		if(args[0].charAt(args[0].length()-1)!=java.io.File.separatorChar){
			args[0]=args[0]+java.io.File.separatorChar;
		}
		if(args[1].equals("")){
			args[1]=args[0]+"MutationsComplete.maf";
		}
		if(args[2].equals("")){
			args[2]=args[0]+"SamplesComplete.txt";
		}
		if(args[3].equals("")){
			args[3]="3";
		}
		if(args[4].equals("")){
			args[4]="1000";
		}
		if(args[5].equals("")){
			args[5]="24";
		}
		if(args[6].equals("")){
			args[6]="100";
		}
		if(args[7].equals("")){
			args[7]="5000";
		}
		if(args[8].equals("")){
			args[8]=args[0]+"src/ComputeDistribution.py";
		}
		if(args[9].equals("")){
			args[9]=args[0]+"Hg19/";
		}
		
		System.out.println("Launching MutPanning with root file "+arguments[0]);
		
		
		System.out.println("Aligning Maf to Hg19 - Step 1");
		AlignHG19.main(new String[]{args[0],args[1],args[2],"0",args[9]});
		System.out.println("Aligning Maf to Hg19 - Step 2");
		AlignHG19.main(new String[]{args[0],args[1],args[2],"1",args[9]});
		System.out.println("Aligning Maf to Hg19 - Step 3");
		AlignHG19.main(new String[]{args[0],args[1],args[2],"2",args[9]});
		System.out.println("Aligning Maf to Hg19 - Step 4");
		AlignHG19.main(new String[]{args[0],args[1],args[2],"3",args[9]});
		
		System.out.println("Determining Count Vectors");
		AffinityCount.main(new String[]{args[0],args[2]});
		AffinityCount_Cosmic.main(new String[]{args[0],args[2]});
		
		System.out.println("Clustering for Each Cancer Type");
		ClusteringEntity.main(new String[]{args[0],args[2],args[5],args[9]});
		
		System.out.println("Clustering PanCancer");
		ClusteringPanCancer.main(new String[]{args[0],args[2],args[3],args[4],args[5],args[9]});
	
		System.out.println("Compute Mutation Rate Clusters");
		ComputeMutationRateClusters.main(new String[]{args[0],args[9]});
	
		System.out.println("Compute Mutation Rate Entities");
		ComputeMutationRateEntities.main(new String[]{args[0],args[2],args[9]});
		
		System.out.println("Count Destructive Mutations");
		CountDestructiveMutations.main(new String[]{args[0],args[1],args[2],args[9]});
		
		System.out.println("Reformat files for CBASE");
		ReformatCBASE.main(new String[]{args[0],args[1],args[2],args[9]});

		
		String[] entities=entities_pancancer(args[2]);
		ArrayList<String>[] samples=samples(args[2],entities);
		int[] count=count(args[0]+"AffinityCounts/AffinityCount.txt",samples);
		
		
		
		
		System.out.println("Execute the parameter estimation for CBASE");
		CBASE_Solutions.main(new String[]{args[0],args[2]});
		
		System.out.println("Compute the mutation count distributions based on CBASE");
		mkdir(args[0]+"CBASE/Distributions/");
		String[] command_step2=new String[entities.length];
		String[] output_files_step2=new String[entities.length];
		for (int i=0;i<entities.length;i++){
			output_files_step2[i]=args[0]+"CBASE/Distributions"+entities[i]+".txt";
			command_step2[i]="python "+args[8]+" \""+args[0]+"CBASE/Counts/CountSilent"+entities[i]+".txt\" \""+args[0]+"CBASE/Counts/Count"+entities[i]+".txt\" \""+args[0]+"CBASE/Parameters_Summary/Parameters"+entities[i]+".txt\" \""+args[0]+"CBASE/Distributions/"+entities[i]+".txt\"";
		}
		execute(command_step2,output_files_step2);
		
		System.out.println("Computing Significance");
		for (int i=0;i<entities.length;i++){
			System.out.println(entities[i]);
			if(samples[i].size()<Integer.parseInt(args[6])||count[i]<Integer.parseInt(args[7])){
				ComputeSignificance_Uniform.main(new String[]{args[0],entities[i],args[9]});
			}
			ComputeSignificance.main(new String[]{args[0],entities[i],args[9]});
		}
		
		System.out.println("Start Filtering of significant genes");
		System.out.println("Filtering Step 1 - Preparing Blat Queries");
	
		
		Filter_Step1.main(new String[]{args[0],args[2],args[9]});
		System.out.println("Performing Blat Queries");
		File[] files_query=new File(args[0]+"PostSignFilter/Queries/").listFiles();
		
		mkdir(args[0]+"PostSignFilter/OutputsBLAT");
		String[] command_blat=new String[files_query.length];
		String[] output_files_blat=new String[files_query.length];
		
		
		for (int i=0;i<files_query.length;i++){			
			output_files_blat[i]=args[0]+"PostSignFilter/OutputsBLAT/Output"+files_query[i].getName().substring(5,files_query[i].getName().length()-3)+".txt";
			command_blat[i]="cd "+args[9]+"SignificanceFilter/"+" && ./blat hg19.2bit "+files_query[i].getAbsolutePath()+" -ooc=11.ooc -out=pslx "+output_files_blat[i];//+" && cd "+args[0];
			
		}
		execute(command_blat,output_files_blat);//,Integer.parseInt(args[5])
		
		System.out.println("Filtering Step 2 - Evaluating Blat Results");
		Filter_Step2.main(new String[]{args[0],args[2]});
		
		System.out.println("Filtering Step 3 - Filtering out Genes");
		Filter_Step3.main(new String[]{args[0],args[2],args[9]});
		
	
	}
	
	
	
	//equally distribute the commands on the no_cpus. make sure that commands are always redistributed as soon as a cpu becomes free
	public static void execute(String[] commands, String[] output_files){
		int no_cpu=24;
		CommandProcess[] processes=new CommandProcess[commands.length];
		for (int i=0;i<processes.length;i++){
			processes[i]=new CommandProcess();
			processes[i].command=commands[i];
			processes[i].output_file=output_files[i];
		}
		
		int no_undone=processes.length;
		int no_running=0;
		while(no_undone>0){
			no_undone=0;
			no_running=0;
			for (int i=0;i<processes.length;i++){
				if(!processes[i].done){
					no_undone++;
				}
				if(processes[i].running){
					no_running++;
				}
			}
			if(no_undone==0){
				break;
			}
			if(no_running<no_cpu){
				for (int i=0;i<processes.length;i++){
					if(!processes[i].done&&!processes[i].running){
						processes[i].start();
						no_running++;
					}
					if(no_running>=no_cpu){
						break;
					}
				}
			}
			try{
				Thread.sleep(3000);
			}
			catch(Exception e){
				System.out.println(e);
			}
		}
	}
	
	
	private static class CommandProcess extends Thread{
		public volatile String command="";
		public volatile String output_file="";
		public volatile boolean done=false;
		public volatile boolean running=false;
		public Process process=null;
		public void run(){
			System.out.println("starting "+command);//+"	"+system_cpu_load()+"	"+process_cpu_load());
			running=true;
			boolean isWindows = System.getProperty("os.name").toLowerCase().startsWith("windows");
			
			if(!new File(output_file).exists()){
				try{
					if (isWindows) {
					    process = Runtime.getRuntime().exec("cmd.exe /c "+command);
					} else {
						String[] commands = { "/bin/bash", "-c", command };
						process = Runtime.getRuntime().exec(commands);
						
					}	
					process.waitFor();
					
				}
				catch(Exception e){
					System.out.println(e);
				}
			}
			
			
			running=false;
			done=true;
			System.out.println("done "+command);//+"	"+system_cpu_load()+"	"+process_cpu_load());
		}
		
		
	}
	
	
	
	public static void mkdir(String file){
		if(!new File(file).exists()){
			new File(file).mkdirs();
		}
	}
	
	public static int[] count(String file,ArrayList<String>[] samples){
		int[] counts=new int[samples.length];
		try{
			FileInputStream in=new FileInputStream(file);
			DataInputStream inn=new DataInputStream(in);
			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
			input.readLine();
			String s="";
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				int count=Integer.parseInt(t[1])+Integer.parseInt(t[2])+Integer.parseInt(t[3])+Integer.parseInt(t[4])+Integer.parseInt(t[5])+Integer.parseInt(t[6]);
				for (int i=0;i<samples.length;i++){
					for (int j=0;j<samples[i].size();j++){
						if(samples[i].get(j).equals(t[0])){
							counts[i]+=count;
						}
					}
				}
			}
			input.close();
		}
		catch(Exception e){
			System.out.println(e);
		}
		return counts;
	}
	
	public static ArrayList<String>[] samples(String file,String[] entities){
		ArrayList<String>[] samples=new ArrayList[entities.length];
		for (int i=0;i<samples.length;i++){
			samples[i]=new ArrayList<String>();
		}
		try{
			FileInputStream in=new FileInputStream(file);
			DataInputStream inn=new DataInputStream(in);
			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
			int[] index_header=index_header(input.readLine().split("	"),index_header_samples);
			String s="";
			while((s=input.readLine())!=null){
				String entity=s.split("	")[index_header[2]];
				String sample=s.split("	")[index_header[1]];
				for (int i=0;i<entities.length;i++){
					if(entities[i].equals("PanCancer")||entities[i].equals(entity)){
						samples[i].add(sample);
					}
				}
			}
			input.close();
		}
		catch(Exception e){
			System.out.println(e);
		}
		return samples;
	}
	
	public static String[] entities(String file){
		try{
			FileInputStream in=new FileInputStream(file);
			DataInputStream inn=new DataInputStream(in);
			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
			int[] index_header=index_header(input.readLine().split("	"),index_header_samples);
			String s="";
			ArrayList<String> aa =new ArrayList<String>();
			while((s=input.readLine())!=null){
				String e=s.split("	")[index_header[2]];
				if(!contains(e,aa)){
					aa.add(e);
				}
			}
			input.close();
			Collections.sort(aa);
			String [] entity=new String[aa.size()];
			for (int i=0;i<aa.size();i++){
				entity[i]=aa.get(i);
			}
			return entity;
		}
		catch(Exception e){
			System.out.println(e);
		}
		return new String[0];
	}
	
	public static String[] entities_pancancer(String file){
		try{
			FileInputStream in=new FileInputStream(file);
			DataInputStream inn=new DataInputStream(in);
			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
			int[] index_header=index_header(input.readLine().split("	"),index_header_samples);
			String s="";
			ArrayList<String> aa =new ArrayList<String>();
			while((s=input.readLine())!=null){
				String e=s.split("	")[index_header[2]];
				if(!contains(e,aa)){
					aa.add(e);
				}
			}
			input.close();
			Collections.sort(aa);
			aa.add("PanCancer");
			String [] entity=new String[aa.size()];
			for (int i=0;i<aa.size();i++){
				entity[i]=aa.get(i);
			}
			return entity;
		}
		catch(Exception e){
			System.out.println(e);
		}
		return new String[0];
	}
	
	public static int[] index_header(String[] header, String[] ideal_header){
		int[] indices=new int[ideal_header.length];
		for (int i=0;i<ideal_header.length;i++){
			int index=-1;
			for (int j=0;j<header.length;j++){
				if(header[j].equals(ideal_header[i])){
					index=j;
					break;
				}
			}
			indices[i]=index;
		}
		return indices;
	}
	
	
	public static boolean contains (String s, ArrayList<String> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).equals(s)){
				return true;
			}
		}
		return false;
	}
}

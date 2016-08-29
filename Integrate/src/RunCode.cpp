/*
 * RunCode.cpp
 *
 *  Created on: Oct 7, 2013
 *      Author: jinzhang
 */

#include "RunCode.h"




RunCode::RunCode() {
	// TODO Auto-generated constructor stub

}

RunCode::~RunCode() {
	// TODO Auto-generated destructor stub
}


void usageBuild()
{
	cout<<endl;
	cout<<"Create directory:"<<endl;
        cout<<endl;
	cout<<"    mkdir directory_to_bwts"<<endl;
	cout<<endl;
	cout<<"Run subcommand mkbwt:"<<endl;
	cout<<endl;
	cout<<"    Integrate mkbwt (options) reference.fasta"<<endl;
	cout<<endl;
	cout<<"    options:"<<endl;
        cout<<endl;
	cout<<"            -mb  integer  :     sequences in the reference fasta that are shorter than this value        default: 10000000"<<endl;
        cout<<"                                are not included in the evaluation of repetitive reads.   "<<endl;
	cout<<"            -dir string   :     directory to store the BWTs.                                             default: ./bwts"<<endl;
	exit(0);
}

int getOptForBuild(int argc, char * argv[], options_t & opt, int & optStart)
{

	for(int i=1;i<argc;i++)
	{
		string tmp=string(argv[i]);
		if(tmp.compare("-mb")==0)
		{
			if(i+1<argc)
			{
				int mb=atoi(argv[i+1]);
				if(mb!=0)
				{
					if(mb>0)
					{
						opt.min_seq_bwt=mb;
					}
					else
					{
						cout<<"Please give a positive integer to -mb"<<endl;
						exit(1);
					}
				}
				else
				{
					string tmp=string(argv[i+1]);
					if(tmp.compare("0")==0)
					{
						opt.min_seq_bwt=0;
					}
					else
					{
						cout<<"Please give an integer after -mb"<<endl;
					}
				}

				if(i+2<argc && i+2>optStart)
					optStart=i+2;
				else
				{
					usageBuild();
				}
			}
			else
			{
				cout<<"Please give an integer after -mb"<<endl;
				exit(1);
			}
		}



		if(tmp.compare("-dir")==0)
		{
			if(i+1<argc)
			{

				//need some code; formal expression

				opt.directoryBWT=argv[i+1];

				if(i+2<argc && i+2>optStart)
					optStart=i+2;
				else
				{
					usageBuild();
				}
			}
			else
			{
				cout<<"Please give a directory after -dir"<<endl;
				exit(1);
			}
		}
	}

	return 0;
}

int RunCode::runBuildBWTs(int argc, char * argv[]) {

	options_t opt;
	opt.min_seq_bwt=10000000;
	opt.directoryBWT="./bwts";

	int opStart=2;

	if(argc<=2)
	{
		usageBuild();
	}

	getOptForBuild(argc,argv,opt,opStart);

	if(opStart<2)
	{
		usageBuild();
	}
	
	if(string(argv[2])=="--help")
	{
		usageBuild();
	}

	float t=clock();
	cout<<"Load ref"<<endl;
	Reference ref;
	ref.setIsInt(0);
	//ref.setIsInt(1);
	ref.test(argv[opStart]);
	cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;


	t=clock();
	cout<<"build BWTs for whole genome"<<endl;
	HitsCounter hc;
	hc.setMinBwtLen(opt.min_seq_bwt);
	hc.getChromBWTs(ref,opt.directoryBWT);

	cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;
	return 0;
}



void usageFusion()
{
    cout<<endl;
    cout<<"Make sure mkbwt has been run."<<endl;
    cout<<endl;
	cout<<"Integrate fusion (options) reference.fasta annotation.txt directory_to_bwt accepted_hits.bam unmapped.bam (dna.tumor.bam dna.normal.bam)\n";
	cout<<endl;
	cout<<"options: -cfn      integer : Cutoff of spanning RNA-Seq reads for fusions with non-canonical"<<endl;
	cout<<"                             exonic boundaries.                                                         default: 3"<<endl;
	cout<<"         -rt       float   : Normal dna / tumor dna ratio. If the ratio is less than"<<endl;
        cout<<"                             this value, then dna reads from the normal dna data set "<<endl;
	cout<<"                             supporting a fusion candidates are ignored.                                default: 0.0"<<endl;
	cout<<"         -minIntra integer : If only having RNA reads, a chimera with two adjacent"<<endl; 
        cout<<"                             genes in order is annotated as intra_chromosomal rather than "<<endl;
	cout<<"                             read_through if the distance of the two genes is longer than"<<endl; 
	cout<<"                             this value.                                                                default: 400000"<<endl; 
	cout<<"         -minW     float   : Mininum weight for the encompassing rna reads on an edge.                  default: 2.0"<<endl; 
	cout<<"         -mb       integer : See subcommand \"mkbwt\"."<<endl; 
	cout<<"                             This value can be larger than used by mkbwt.                               default: 10000000"<<endl;
	cout<<"         -minDel   int     : minimum size of a deletion that can cause a fusion.                        default: 5000"<<endl;
	cout<<"         -reads    string  : File to store all the reads.                                               default: reads.txt"<<endl;
	cout<<"         -sum      string  : File to store summary.                                                     default: summary.tsv"<<endl;
	cout<<"         -ex       string  : File to store exons for fusions with canonical exonic boundaries.          default: exons.tsv"<<endl;
	cout<<"         -bk       string  : File to store breakpoints                                                  default: breakpoints.tsv"<<endl;
	cout<<"         -vcf      string  : File to store breakpoints in vcf format                                    default: bk_sv.vcf"<<endl;

        cout<<"         -bedpe   string  : File to store all fusions in SMC-RNA bedpe format                           default: fusions.bedpe"<<endl;   
	cout<<"         -bacc     integer : max difference between spanning reads and annotation to decide canonical.  default: 1"<<endl;
        cout<<"         -largeNum integer : if a gene shows greater or equal to this number, remove it from results.   default: 4"<<endl;
        cout<<"         -sample   string  : sample name                                                                default: sample"<<endl;
	cout<<endl;
	cout<<"This version of Integrate works in the following situations:"<<endl;
	cout<<"(1)having rna tumor, dna tumor, dna normal"<<endl;
	cout<<"(2)having rna tumor, dna tumor"<<endl;
	cout<<"(3)having rna tumor"<<endl;
	//cout<<"dna bam files should be followed by its library insert size and standard deviation. If RG lines are provided, then insert size is ignored, but a value still has to be provided."<<endl;
	cout<<endl;
	cout<<"Integrate will only use sequences in reference.fasta. "<<endl;
	cout<<"Chr names with and without \"chr\" are regarded as the same, e.g. chr1 = 1."<<endl;
	cout<<"The rna and dna bams can be from alignments mapped to different reference files with different order of the sequences and their names with or without \"chr\". However, The versions should be the same, e.g. hg19. (Also, the same as in annotation.)"<<endl;
	cout<<"The tumor and normal dna bams should be mapped to the same reference file."<<endl;
	cout<<endl;
	cout<<"For rna tumor: accepted_hits.bam is a bam file containing mapped rna reads. unmapped.bam is a bam contains the not mapped rna reads. If they have been merged into one bam, just use merged.bam twice in the command line."<<endl;
	cout<<endl;
	cout<<"For dna bams: If solt-clips are provided, then Integrate is trying to search rearrangement breakpoints, otherwise, only paired reads may be included in the analysis."<<endl;
	cout<<endl;
	cout<<"If having rna normal only or having both rna and dna normal data sets. These data sets can be run to find non somatic events."<<endl;
	cout<<"e.g. Integrate fusion -normal (options) reference.fasta annotation.txt directory_to_bwt accepted_hits.normal.bam unmapped.normal.bam (dna.normal.bam)"<<endl;
	exit(0);
}

int getOptForFusion(int argc, char * argv[], options_t & opt, int & opStart)
{

	for(int i=1;i<argc;i++)
	{




		string tmp=string(argv[i]);

		if(tmp.compare("-mb")==0)
		{
			if(i+1<argc)
			{
				int mb=atoi(argv[i+1]);
				if(mb!=0)
				{
					if(mb>0)
					{
						opt.min_seq_bwt=mb;
					}
					else
					{
						cout<<"Please give a positive integer to -mb"<<endl;
						exit(1);
					}
				}
				else
				{
					string tmp=string(argv[i+1]);
					if(tmp.compare("0")==0)
					{
						opt.min_seq_bwt=0;
					}
					else
					{
						cout<<"Please give an integer after -mb"<<endl;
					}
				}

				if(i+2<argc && i+2>opStart)
					opStart=i+2;
				else
				{
					usageBuild();
				}
			}
			else
			{
				cout<<"Please give an integer after -mb"<<endl;
				exit(1);
			}
		}


		if(tmp.compare("-cfn")==0)
		{
			if(i+1<argc)
			{
				int mb=atoi(argv[i+1]);


				if(mb>0)
				{
					opt.cfn=mb;
				}
				else
				{
					cout<<"Please give an integer >= 1 to -cfn"<<endl;
					exit(1);
				}


				if(i+2<argc && i+2>opStart)
					opStart=i+2;
				else
				{
					usageFusion();
				}
			}
			else
			{
				cout<<"Please give an integer after -cfn"<<endl;
				exit(1);
			}
		}

                if(tmp.compare("-minDel")==0)
                {
                        if(i+1<argc)
                        {
                                int mb=atoi(argv[i+1]);


                                if(mb>0)
                                {
                                        opt.minDel=mb;
                                }
                                else
                                {
                                        cout<<"Please give an integer >= 1 to -minDel"<<endl;
                                        exit(1);
                                }


                                if(i+2<argc && i+2>opStart)
                                        opStart=i+2;
                                else
                                {
                                        usageFusion();
                                }
                        }
                        else
                        {
                                cout<<"Please give an integer after -minDel"<<endl;
                                exit(1);
                        }
                }


		if(tmp.compare("-rt")==0)
		{
			if(i+1<argc)
			{
				double mb=atof(argv[i+1]);

				if(mb>=0.0 && mb<=1.0 )
				{
					opt.rt=mb;
				}
				else
				{
					cout<<"Please give a ration >=0.0 and <=1.0 to -rt"<<endl;
					exit(1);
				}


				if(i+2<argc && i+2>opStart)
					opStart=i+2;
				else
				{
					usageFusion();
				}
			}
			else
			{
				cout<<"Please give an ratio after -rt"<<endl;
				exit(1);
			}
		}


		if(tmp.compare("-minW")==0)
		{
			if(i+1<argc)
			{
				double mb=atof(argv[i+1]);

				if(mb<1.0 )
				{
					cout<<"Please give a value >=1.0 to -minW"<<endl;
					exit(1);
				}
				else
				{
					opt.minW=mb;
				}


				if(i+2<argc && i+2>opStart)
					opStart=i+2;
				else
				{
					usageFusion();
				}
			}
			else
			{
				cout<<"Please give a value >=1.0 after -minW"<<endl;
				exit(1);
			}
		}


		if(tmp.compare("-minIntra")==0)
		{
			if(i+1<argc)
			{
				int mb=atoi(argv[i+1]);

				if(mb<=0)
				{
					cout<<"Please give a positive integer to -minIntra"<<endl;
					exit(1);
				}
				else if(mb<=70000)
				{
					cout<<"Warning: the value of -minIntra is "<<mb<<". This value may be very small, rna only events may be annotated as intra_chromosomal"<<endl;
					opt.rt=mb;
				}
				else
				{
					opt.minIntra=mb;
				}


				if(i+2<argc && i+2>opStart)
					opStart=i+2;
				else
				{
					usageFusion();
				}
			}
			else
			{
				cout<<"Please give an Integer after -minIntra"<<endl;
				exit(1);
			}

		}


		if(tmp.compare("-reads")==0)
		{
			if(i+1<argc)
			{
				//need some code; formal expression

				opt.fileRead=argv[i+1];

				if(i+2<argc && i+2>opStart)
					opStart=i+2;
				else
				{
					usageFusion();
				}
			}
			else
			{
				cout<<"Please give a file name after -reads"<<endl;
				exit(1);
			}
		}


		if(tmp.compare("-sum")==0)
		{
			if(i+1<argc)
			{

				//need some code; formal expression

				opt.fileSum=argv[i+1];

				if(i+2<argc && i+2>opStart)
					opStart=i+2;
				else
				{
					usageFusion();
				}
			}
			else
			{
				cout<<"Please give a file name after -sum"<<endl;
				exit(1);
			}
		}

		if(tmp.compare("-ex")==0)
		{
			if(i+1<argc)
			{

				//need some code; formal expression

				opt.fileEx=argv[i+1];

				if(i+2<argc && i+2>opStart)
					opStart=i+2;
				else
				{
					usageFusion();
				}
			}
			else
			{
				cout<<"Please give a file name after -ex"<<endl;
				exit(1);
			}
		}

        
        if(tmp.compare("-bk")==0)
		{
			if(i+1<argc)
			{
                
				//need some code; formal expression
                
				opt.bkFile=argv[i+1];
                
				if(i+2<argc && i+2>opStart)
					opStart=i+2;
				else
				{
					usageFusion();
				}
			}
			else
			{
				cout<<"Please give a file name after -bk"<<endl;
				exit(1);
			}
		}
       
/*        if(tmp.compare("-bedpe")==0)
                {
                        if(i+1<argc)
                        {

                                //need some code; formal expression

                                opt.bkFileBEDPE=argv[i+1];

                                if(i+2<argc && i+2>opStart)
                                        opStart=i+2;
                                else
                                {
                                        usageFusion();
                                }
                        }
                        else
                        {
                                cout<<"Please give a file name after -bedpe"<<endl;
                                exit(1);
                        }
                }

*/
        if(tmp.compare("-vcf")==0)
                {
                        if(i+1<argc)
                        {

                                //need some code; formal expression

                                opt.bkFileVCF=argv[i+1];

                                if(i+2<argc && i+2>opStart)
                                        opStart=i+2;
                                else
                                {
                                        usageFusion();
                                }
                        }
                        else
                        {
                                cout<<"Please give a file name after -vcf"<<endl;
                                exit(1);
                        }
                }

        if(tmp.compare("-junc")==0)
                {
                        if(i+1<argc)
                        {

                                //need some code; formal expression

                                opt.fileJunction=argv[i+1];

                                if(i+2<argc && i+2>opStart)
                                        opStart=i+2;
                                else
                                {
                                        usageFusion();
                                }
                        }
                        else
                        {
                                cout<<"Please give a file name after -junc"<<endl;
                                exit(1);
                        }
                }

        if(tmp.compare("-peptide")==0)
                {
                        if(i+1<argc)
                        {

                                //need some code; formal expression

                                opt.filePeptide=argv[i+1];

                                if(i+2<argc && i+2>opStart)
                                        opStart=i+2;
                                else
                                {
                                        usageFusion();
                                }
                        }
                        else
                        {
                                cout<<"Please give a file name after -peptide"<<endl;
                                exit(1);
                        }
                }

        if(tmp.compare("-bedpe")==0)
                {
                        if(i+1<argc)
                        {

                                //need some code; formal expression

                                opt.fileSmcRna=argv[i+1];

                                if(i+2<argc && i+2>opStart)
                                        opStart=i+2;
                                else
                                {
                                        usageFusion();
                                }
                        }
                        else
                        {
                                cout<<"Please give a file name after -smcRna"<<endl;
                                exit(1);
                        }
                }

        if(tmp.compare("-sample")==0)
                {
                        if(i+1<argc)
                        {

                                //need some code; formal expression

                                opt.sample_name=argv[i+1];

                                if(i+2<argc && i+2>opStart)
                                        opStart=i+2;
                                else
                                {
                                        usageFusion();
                                }
                        }
                        else
                        {
                                cout<<"Please give a sample name after -sample"<<endl;
                                exit(1);
                        }
                }

		if(tmp.compare("-normal")==0)
		{
			opt.isRunningNormal=1;
			opStart=i+1;
		}


		if(tmp.compare("-bacc")==0)
                {
                        if(i+1<argc)
                        {
                                int mb=atoi(argv[i+1]);


                                if(mb>=0)
                                {
                                        opt.bacc=mb;
                                }
                                else
                                {
                                        cout<<"Please give an integer >= 0 to -bacc"<<endl;
                                        exit(1);
                                }


                                if(i+2<argc && i+2>opStart)
                                        opStart=i+2;
                                else
                                {
                                        usageFusion();
                                }
                        }
                        else
                        {
                                cout<<"Please give an integer after -bacc"<<endl;
                                exit(1);
                        }
                }
		
		if(tmp.compare("-largeNum")==0)
                {
                        if(i+1<argc)
                        {
                                int mb=atoi(argv[i+1]);


                                if(mb>=2)
                                {
                                        opt.largeNum=mb;
                                }
                                else
                                {
                                        cout<<"Please give an integer >= 2 to -largeNum"<<endl;
                                        exit(1);
                                }


                                if(i+2<argc && i+2>opStart)
                                        opStart=i+2;
                                else
                                {
                                        usageFusion();
                                }
                        }
                        else
                        {
                                cout<<"Please give an integer after -largeNum"<<endl;
                                exit(1);
                        }
                }


	}
	return 0;
}


int RunCode::runFindFusions(int argc, char * argv[]) {




		options_t opt;
		int opStart=2;
		
		opt.cfn=3;
		opt.rt=0.0;
		opt.minIntra=400000;
		opt.minW=2.0;
		opt.min_seq_bwt=10000000;
		opt.isRunningNormal=0;
		opt.bacc=1;
		opt.largeNum=4;
                opt.minDel=5000;

		if(argc<=2)
		{
			usageFusion();
		}

		//getOptForFusion(argc,argv,opt,opStart);

		if(opStart<2)
		{
			usageFusion();
		}

		if(opt.isRunningNormal==1)
		{
			opt.fileRead="reads_normal.tsv";
			opt.fileSum="summary_normal.tsv";
			opt.fileEx="exons_normal.tsv";
            		opt.bkFile="breakpoints_normal.tsv";
			opt.bkFileBEDPE="bk_fusion_normal.bedpe";
			opt.bkFileVCF="bk_sv_normal.vcf";
                        opt.fileJunction="junctions_normal.bedpe";
                        opt.filePeptide="peptides_normal.bedpe";
                        opt.fileSmcRna="fusions_normal..bedpe";
 			opt.sample_name="sample_normal";
		}
		else
		{
			opt.fileRead="reads.txt";
			opt.fileSum="summary.tsv";
			opt.fileEx="exons.tsv";
            		opt.bkFile="breakpoints.tsv";
			opt.bkFileBEDPE="bk_fusion.bedpe";//depleted
			opt.bkFileVCF="bk_sv.vcf";
                        opt.fileJunction="junctions.bedpe";
                        opt.filePeptide="peptides.bedpe";
                        opt.fileSmcRna="fusions.bedpe";	
			opt.sample_name="sample_tumor";
 		}

		getOptForFusion(argc,argv,opt,opStart);

		Result result;


		int insert1,std1;
		int insert2,std2;




		result.setIndi(2);
		if(argc-opStart==7)
		{
			/*

			insert1=atoi(argv[opStart+6]);
			std1=atoi(argv[opStart+7]);

			if(insert1<=0 || std1<=0)
			{
				usageFusion();
			}


			insert2=atoi(argv[opStart+9]);
			std2=atoi(argv[opStart+10]);

			if(insert2<=0 || std2<=0)
			{
				usageFusion();
			}

			*/

			result.setIndi(3);

		}else if(argc-opStart==5)
		{
			result.setIndi(1);
		}
		else if(argc-opStart==6)
		{
			/*
			insert1=atoi(argv[opStart+6]);
			std1=atoi(argv[opStart+7]);

			if(insert1<=0 || std1<=0)
			{
				usageFusion();
			}
			*/
		}
		else
		{
			usageFusion();
		}




		myFind2 mf2;
		mf2.setMaxdiff(2);
	    mf2.setStateLen(1000000);
	    mf2.create();

		float t=clock();
		cout<<"Loading reference..."<<endl;
		Reference ref;
		ref.setIsInt(0);
		ref.test(argv[opStart]);
		cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;


		TidHandler th;
                th.setRefTid(ref);


                Gene g;

                t=clock();
                cout<<"Loading genes..."<<endl;
                g.loadGenesFromFile(argv[opStart+1],th);
                g.setGene();
                g.allocate();
                cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;



		t=clock();
		cout<<"Loading BWTs for chromosomes..."<<endl;
		HitsCounter hc;
		hc.setMinBwtLen(opt.min_seq_bwt);
		hc.loadChromBWTs(ref,argv[opStart+2]);
		cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;


		cout<<"Handling RNA data...\n"<<endl;

		MyBamHeader rnabh;
		rnabh.run2(argv[opStart+3]);//300 and 100 is not used
		th.setRNAAndRef(rnabh);


		t=clock();
		cout<<"\nGetting graph by encompassing RNA reads..."<<endl;
		Rna rna;
		rna.getGraph(argv[opStart+3],th, g, hc);
		//cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;


		//t=clock();
		//cout<<"combine records"<<endl;
		rna.cbEncompassRcs(g);
		cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;


		t=clock();
		cout<<"Reducing Graph by removing encompassing pair with a hit in one gene mapped in BAM ..."<<endl;
		rna.reduceGraph(argv[opStart+3],g,th);
		cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;


		MyBamWrap mbw;
		mbw.mySamOpen(argv[opStart+3]);
		mbw.myGetIndex(argv[opStart+3]);


		t=clock();
		cout<<"Getting BWTs for genes in the graph..."<<endl;
		rna.runGetGeneBWT(g,ref);
		cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;


		t=clock();
		cout<<"Adding spanning RNA reads from BAM file..."<<endl;
		rna.readTopHat2(argv[opStart+3],th, g,hc, ref,mf2);
		rna.handleTmpTopHatSplits(hc,g);
		cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;

		cout<<"Handle secondary split reads..."<<endl;
		rna.readSTAR(argv[opStart+3],th, g,hc, ref,mf2);
		cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;

		t=clock();
		cout<<"Computing number of reads on edges..."<<endl;
		rna.computeWeights1(g);
		cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;

		t=clock();
		cout<<"Removing edges with too few reads..."<<endl;
		rna.reduceGraph2(g,opt.minW);
		cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;


		t=clock();
		cout<<"Removing edges caused by repetitive encompassing reads by mapping them to both genes... "<<endl;
		rna.reduceInOneGene(g,mf2);
		rna.reduceGraphEmptyEn(g);
		cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;

		t=clock();
		cout<<"Estimating how repetitive of the remaining reads..."<<endl;
		rna.computeNumCopy(hc);
		cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;


		t=clock();
		cout<<"Computing weights for edges..."<<endl;
		rna.computeWeights1(g);
		cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;

		t=clock();
		cout<<"Removing edges with weight less than "<<opt.minW<<" ..."<<endl;
		rna.reduceGraph2(g,opt.minW);
		cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;

		cout<<"Getting hard cliped secondary reads' sequences ..."<<endl;
		rna.getHardClipReads(g,mbw,th);
		cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;

		t=clock();
		cout<<"Mapping soft or hard clipped reads and reads that could be partially correctly mapped from BAM"<<endl;
		rna.traversePartialRight(g, mbw, th, hc, mf2);
		rna.handleTmpTopHatSplits(hc,g);
		cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;

		t=clock();
		cout<<"Getting anchors from unaligned reads"<<endl;
		rna.getAnchors(g,mbw,th, hc);
		cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;

		t=clock();
		cout<<"Mapping the unmapped reads as split reads..."<<endl;
		rna.mapPartialSplitBWT(argv[opStart+4],th,g,hc,mf2);
		cout<<"Total time : "<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;


		t=clock();
		cout<<"Handling mapped split reads..."<<endl;
		rna.handleSpHits();
		//cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;


		//t=clock();
		//cout<<"cluster "<<endl;
		rna.traverseCluster(g,opt.cfn,opt.bacc);
		//cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;



		//t=clock();
		//cout<<"print "<<endl;
		rna.traversePrint(g,ref,result,opt.minIntra,opt.bacc);
		cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;

		result.checkALLPrime();

		if(result.getIndi()>1)
		{

		t=clock();
		cout<<"Processing DNA Tumor...\n"<<endl;
		MyBamHeader dnabh;
		dnabh.run2(argv[opStart+5]);
		th.setDNAAndRef(dnabh);
		//cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;

		//t=clock();
		//cout<<"run"<<endl;
		Dna dna;
		dna.onlyDNAByResult(argv[opStart+5],g,th,dnabh,ref,result,opt.minDel, 0);
		cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;

			if(result.getIndi()>2)
			{
				t=clock();
				cout<<"Processing DNA Normal...\n"<<endl;
				MyBamHeader dnabh2;
				dnabh2.run2(argv[opStart+6]);
				//th.setDNAAndRef(dnabh);// current it is same?
				//cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;

				//t=clock();
				//cout<<"run"<<endl;
				Dna dna2;
				dna2.onlyDNAByResult(argv[opStart+6],g,th,dnabh2,ref,result, opt.minDel, 1);
				cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;
			}

		}
		t=clock();
		cout<<"Printing results..."<<endl;
		result.getTiers(opt.rt);
		result.printSummary(opt.fileSum, g, opt.isRunningNormal, opt.largeNum);
		result.printAllResult(opt.fileRead,ref,opt.isRunningNormal);
		result.printExons(opt.fileEx,g,ref,opt.isRunningNormal, opt.bkFile, opt.bkFileBEDPE, opt.bkFileVCF, argv[opStart], opt.sample_name);
		//Dec 7, 2015
                result.getAllJunctionsStep4(g,ref);
                result.getAllJunctionsStep5(g,ref);
                result.getAllJunctionsStep6(opt.fileSmcRna,g,ref);
		cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;

		cout<<"Done. Please refer to summary in "<<opt.fileSum<<"."<<endl;
                cout<<"      Encompassing and Spanning reads are in "<<opt.fileRead<<"."<<endl;
		cout<<"      To design probes for validation, please refer to "<<opt.fileEx<<"."<<endl;
		cout<<"      For beakpoints, please refer to "<<opt.bkFile<<", "<<opt.fileSmcRna<<", and "<<opt.bkFileVCF<<"."<<endl;
                cout<<"      Note:"<<opt.fileSmcRna<<" is in SMC-RNA bedpe format."<<endl; 
		return 0;
}



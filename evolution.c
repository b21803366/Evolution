#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct Gene																				
{
	int data;
	struct Gene* link;
};

struct Chromesome
{
	struct Gene* genes;
	struct Chromesome* link;
	int fitness;
	double rank;
};

void assignRank(struct Chromesome** populationP,int POP_SIZE);									//assign ranks
void assignFitness(struct Chromesome** populationP,int POP_SIZE,int MAX_GEN);					//assign fitnesses
int rankToIndex(struct Chromesome* population,int popHead,int POP_SIZE,int rank);  				//convert link list rank to array index
void sortLinkedList(struct Chromesome** populationP,int POP_SIZE,int MAX_GEN,int* popHeadP);	//sorts members with links
void printPopulation(struct Chromesome* population,int popHead);								//prints population
int calFitness(struct Chromesome chromesome,int MAX_GEN);										//calculates fitness of a chromesome
char** split(char* string,char spliter);														//splits string according to spliter
int len(char* string);																			//length of a string (this function for dev.cs)
char* readFile(char* fileName);																	//reads all of the file
int fileLen(char* fileName);																	//returns char count of a file

int main(int argc,char *argv[])
{
	
	if(argc!=4)																					//argument count check
	{
		printf("Not enough arguments!\n");
		return;
	}
	const int PROB_SIZE = atoi(argv[1]);														//assigning constants
	const int POP_SIZE = atoi(argv[2]);
	const int MAX_GEN = atoi(argv[3]);
	
	/*
	const int PROB_SIZE = 10;
	const int POP_SIZE = 8;
	const int MAX_GEN = 10;
	*/
	int i,j,k,temp,sel1,sel2,mut,xover1,xover2,popHead;
	
	struct Chromesome* population = (struct Chromesome*)malloc(POP_SIZE*sizeof(struct Chromesome));	//creating population according to POP_SIZE
	
	struct Chromesome bestChromesome;
	
	//freopen("output.txt","w",stdout);
	
	for(i=0;i<(POP_SIZE-1);i++)																		//default linking for chromesomes
	{
		population[i].link=(&population[i+1]);
	}
	population[(POP_SIZE-1)].link = NULL;
	popHead = 0;
	
	for(i=0;i<POP_SIZE;i++)
	{
		population[i].genes = (struct Gene*)malloc((MAX_GEN)*sizeof(struct Gene));					//creating genes according to MAX_GEN
		
		for(j=0;j<(MAX_GEN-1);j++)																	//default linking for genes
		{
			population[i].genes[j].link=(&population[i].genes[j+1]);			
		}
		population[i].genes[MAX_GEN-1].link = NULL;
		
		for(j=0;j<MAX_GEN;j++)
		{
			population[i].genes[j].data = atoi(split(split(readFile("population"),'\n')[i],':')[j]);	//assigning genes datas
		}
	}
	
	printf("GENERATION: 0\n");
	
	assignFitness(&population,POP_SIZE,MAX_GEN);												//after assigning datas, calculating fitnesses and assigning them
	assignRank(&population,POP_SIZE);															//calculating and assigning ranks
	sortLinkedList(&population,POP_SIZE,MAX_GEN,&popHead);										//sort chromesomes according to rank values
	
	printPopulation(population,popHead);														//prints population
	
	bestChromesome = population[popHead];														//assigning best chromesome
	
	printf("Best chromosome found so far: ");													//printing best chromesome
	printf("%d",bestChromesome.genes[0].data);
	for(j=1;j<MAX_GEN;j++)
	{
		printf(":%d",bestChromesome.genes[j].data);
	}
	printf(" -> %d\n",bestChromesome.fitness);
	
	for(i=0;i<PROB_SIZE;i++)
	{
		mut = atoi(split(readFile("mutate"),'\n')[i]);											//takes mutate value in i. line
		
		xover1 = atoi(split(split(readFile("xover"),'\n')[i],':')[0]);							//takes xover values in i. line
		xover2 = atoi(split(split(readFile("xover"),'\n')[i],':')[1]);
		
		for(j=0;j<(POP_SIZE/2);j++)
		{
			sel1 = atoi(split(split(split(readFile("selection"),'\n')[i],' ')[j],':')[0]);		//takes POP_SIZE/2 different selection in i. line
			sel1 = rankToIndex(population,popHead,POP_SIZE,sel1-1);
			
			sel2 = atoi(split(split(split(readFile("selection"),'\n')[i],' ')[j],':')[1]);
			sel2 = rankToIndex(population,popHead,POP_SIZE,sel2-1);
			
			for(k=(xover1-1);k<=(xover2-1);k++)													//xovers
			{
				temp = population[sel1].genes[k].data;
				population[sel1].genes[k].data = population[sel2].genes[k].data;
				population[sel2].genes[k].data = temp;
			}
			
			population[sel1].genes[mut-1].data = (1-population[sel1].genes[mut-1].data);		//mutations
			population[sel2].genes[mut-1].data = (1-population[sel2].genes[mut-1].data);			
			
		}
		printf("GENERATION: %d\n",(i+1));
		
		assignFitness(&population,POP_SIZE,MAX_GEN);											//calculating and assigning fitnesses
		assignRank(&population,POP_SIZE);														//calculating and assigning ranks
		sortLinkedList(&population,POP_SIZE,MAX_GEN,&popHead);									//sort chromesomes according to rank values
		
		printPopulation(population,popHead);													//prints population
		if(bestChromesome.fitness > population[popHead].fitness)
			bestChromesome = population[popHead];
		
		printf("Best chromosome found so far: ");
		
		printf("%d",bestChromesome.genes[0].data);
		for(j=1;j<MAX_GEN;j++)
		{
			printf(":%d",bestChromesome.genes[j].data);
		}
		printf(" -> %d\n",bestChromesome.fitness);
	}
	
	return 0;
}

void assignFitness(struct Chromesome** populationP,int POP_SIZE,int MAX_GEN)
{
	int i;
	for(i=0;i<POP_SIZE;i++)
	{
		(*populationP)[i].fitness = calFitness((*populationP)[i],MAX_GEN);
	}
}

void assignRank(struct Chromesome** populationP,int POP_SIZE)
{
	int i;
	double bestFitness = 1024;
	double worstFitness = 0;
	
	for(i=0;i<POP_SIZE;i++)
	{
		if((*populationP)[i].fitness < bestFitness)											//finds best fitness value
		{
			bestFitness = (*populationP)[i].fitness;
		}
		if((*populationP)[i].fitness > worstFitness)										//finds worst fitness value
		{
			worstFitness = (*populationP)[i].fitness;
		}
	}
	
	if(bestFitness == worstFitness)
	{
		for(i=0;i<POP_SIZE;i++)																//default
		{
			(*populationP)[i].rank = -1;
		}
	}
	else
	{
		for(i=0;i<POP_SIZE;i++)
		{
			(*populationP)[i].rank = ((double)((*populationP)[i].fitness)-bestFitness)/(worstFitness-bestFitness);
		}
	}
}

void sortLinkedList(struct Chromesome** populationP,int POP_SIZE,int MAX_GEN,int* popHeadP)
{
	int i,j,k,preIndex;
	int minIndex = 0;
	double minValue = 1;
	int* assigned = (int*)malloc(POP_SIZE*sizeof(int));
	int assignCount=0;
	int checked;
	
	for(i=0;i<POP_SIZE;i++)
	{
		assigned[i] = (-1);
	}
	
	for(i=0;i<POP_SIZE;i++)															//base case
	{
		if((*populationP)[i].rank < minValue)
		{
			minIndex = i;
			minValue = (*populationP)[i].rank;
		}
	}
	
	(*popHeadP) = minIndex;
	preIndex = minIndex;
	assignCount++;
	assigned[assignCount-1] = minIndex;
	
	for(i=0;i<(POP_SIZE-1);i++)																		//link according to rank
	{
		minValue = 1024;
		for(j=0;j<POP_SIZE;j++)
		{
			checked = 0;
			for (k=0;k<assignCount;k++)
			{
				if(j == assigned[k])
				{
					checked = 1;
				}
			}
			if( (checked == 1) && ((*populationP)[j].rank == (*populationP)[preIndex].rank) )		//check if already sorted
			{
				continue;
			}
			else if( ((*populationP)[j].rank < minValue) && ((*populationP)[j].rank >= (*populationP)[preIndex].rank) )
			{
				minIndex = j;
				minValue = (*populationP)[minIndex].rank;
			}
		}
		(*populationP)[preIndex].link = (&(*populationP)[minIndex]);
		preIndex = minIndex;
		assignCount++;
		assigned[assignCount-1] = minIndex;
	}
	(*populationP)[minIndex].link = NULL;
	
}

void printPopulation(struct Chromesome* population,int popHead)
{
	struct Chromesome currentChromesome = population[popHead];
	struct Gene currentGene;
	while(1)
	{
		currentGene = currentChromesome.genes[0];
		while(1)
		{
			printf("%d ",currentGene.data);
			if(currentGene.link == NULL)
				break;
			else
				currentGene = (*currentGene.link);
		}
		printf("-> %d",currentChromesome.fitness);
		printf("\n");
		if(currentChromesome.link == NULL)
			break;
		else
			currentChromesome = (*currentChromesome.link);
	}
}

int rankToIndex(struct Chromesome* population,int popHead,int POP_SIZE,int rank)
{
	int i,index;
	if(rank==0)
	{
		return popHead;	
	}
		
	
	struct Chromesome currentChrome = population[popHead];
	for(i=0;i<rank;i++)
	{
		currentChrome = (*currentChrome.link);
	}
	
	for(i=0;i<POP_SIZE;i++)
	{
		if(population[i].fitness == currentChrome.fitness)
			return i;
	}
}

int calFitness(struct Chromesome chromesome,int MAX_GEN)
{
	int i;
	int fitness=0;
	
	for(i=0;i<MAX_GEN;i++)
	{
		fitness += (pow(2,(MAX_GEN-i-1))*chromesome.genes[i].data);
	}
	return fitness;
}

char** split(char* string,char spliter)
{
	int i,charCount=0;
	int partCount=0;
	char** parts = (char**)malloc(1 * sizeof(char*));
	parts[0] = (char*)malloc(1 * sizeof(char)); 
	
	for(i=0;i<len(string);i++)
	{
		if(string[i]==spliter)
		{
			if(string[i+1]==spliter)
				continue;
			else if(string[i+1]=='\0')
				break;
				
			charCount++;
			if(charCount!=1)
				parts[partCount] = (char*)realloc(parts[partCount],charCount*sizeof(char));
			parts[partCount][charCount-1] = '\0';
			
			partCount++;
			charCount=0;
			
			parts = (char**)realloc(parts,(partCount+1)*sizeof(char*));
			parts[partCount] = (char*)malloc(1 * sizeof(char));
		}
		else
		{
			charCount++;
			if(charCount!=1)
				parts[partCount] = (char*)realloc(parts[partCount],charCount*sizeof(char));
			parts[partCount][charCount-1] = string[i];
		}
	}
	charCount++;
	if(charCount!=1)
		parts[partCount] = (char*)realloc(parts[partCount],charCount*sizeof(char));
	parts[partCount][charCount-1] = '\0';
	
	partCount++;
	return parts;
}

int len(char* string)
{
	int i;
	for(i=0;;i++)
	{
		if(string[i]=='\0')
			break;
	}
	return i;
}

char* readFile(char* fileName)
{
	char ch;
	int charCount = fileLen(fileName);
	char* fileStr = (char*)malloc((charCount+1)*sizeof(char));
	FILE* file = fopen(fileName,"r");
	int i;
	
	for(i=0;i<charCount;i++)
	{
		ch = fgetc(file);
		if(ch==EOF)
		{
			break;
		}
		fileStr[i]=ch;
	}
	fileStr[charCount]='\0';
	
	fclose(file);
	return fileStr;
}

int fileLen(char* fileName)
{
	char ch;
	FILE* file = fopen(fileName,"r");
	int charCount=0;
	
	while(1)
	{
		ch = fgetc(file);
		if(ch==EOF)
		{
			break;
		}
		charCount++;
	}
	fclose(file);
	return charCount;
}

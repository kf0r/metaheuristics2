#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <stack>
#include <stack>
#include <queue>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <random>
#include <chrono>

using namespace std;

random_device r;
mt19937 g(r());

struct indexPair{
    int first;
    int second;
};

struct Node{
    int x;
    int y;
    int id;
    bool visited;
};

struct TreeNode{
    vector<TreeNode*> neighbours;
    int index;
    bool visited;
};

class Tree{
    TreeNode** nodes;

public:
    int size;
    Tree(int size){
        nodes = new TreeNode*[size];
        for(int i =0; i<size; i++){
            nodes[i] = new TreeNode;
        }
        this->size=size;
    }

    ~Tree(){
        for (int i=0; i<size; i++){
            delete nodes[i];
        }
        delete[] nodes;
    }

    void addNode(int n1, int n2){
        nodes[n1]->neighbours.push_back(nodes[n2]);
        nodes[n2]->neighbours.push_back(nodes[n1]);
        nodes[n1]->visited=false;
        nodes[n2]->visited=false;
        nodes[n1]->index=n1;
        nodes[n2]->index=n2;
    }

    void resetTree(){
        stack<int> stack;
        stack.push(0);

        while(!stack.empty()){
            int s = stack.top();
            stack.pop();
            if(nodes[s]->visited){
                nodes[s]->visited=false;
            }

            for(int i = 0; i<nodes[s]->neighbours.size(); i++){
                if(nodes[s]->neighbours[i]->visited){
                    stack.push(nodes[s]->neighbours[i]->index);
                }
            }
        }
    }

    vector<int> dfs(int start){
        vector<int> cycle;
        stack<int> stack;
        stack.push(start);

        while(!stack.empty()){
            int s = stack.top();
            stack.pop();
            if(!nodes[s]->visited){
                nodes[s]->visited=true;
                cycle.push_back(s);
            }

            for(int i = 0; i<nodes[s]->neighbours.size(); i++){
                if(!nodes[s]->neighbours[i]->visited){
                    stack.push(nodes[s]->neighbours[i]->index);
                }
            }
        }
        resetTree();
        return cycle;
    }
};

class Graph{
    Node* nodes;
    

    int calcDist(Node n1, Node n2){
        int dx = n1.x-n2.x;
        int dy = n1.y-n2.y;
        int dist2 = dx*dx+dy*dy;
        return (int)round(sqrt(dist2));
    }

    bool isValidEdge(int u, int v){
        if (u==v){
            return false;
        }
        if(nodes[u].visited&&nodes[v].visited){
            return false;
        }
        if(nodes[u].visited==false&&nodes[v].visited==false){
            return false;
        }
        return true;
    }

public:
    int size;
    int** adjacencyMatrix;
    Tree* tree;
    Graph(int size){
        this->size = size;
        adjacencyMatrix = new int*[size];
        for (int i = 0; i < size; ++i) {
            adjacencyMatrix[i] = new int[size];
        }
        nodes = new Node[size];
        tree = new Tree(size);
    }

    void loadNode(int x, int y, int id){
        nodes[id].id = id;
        nodes[id].x = x;
        nodes[id].y = y;
        nodes[id].visited=false;
    }

    void prepareGraph(){
        for(int i=0; i<size;i++){
            for(int j=0; j<size;j++){
                adjacencyMatrix[i][j]=calcDist(nodes[i], nodes[j]);
            }
        }
    }

    ~Graph() {
        for (int i = 0; i < size; ++i) {
            delete[] adjacencyMatrix[i];
        }
        delete[] adjacencyMatrix;
    }

    int monteCarloCycle(int repeats){
        int* cycle = new int[size];
        for (int i = 0; i<size; i++){
            cycle[i] = i;
        }

        int cost;
        int min = INT32_MAX;
        for(int j=0; j<repeats; j++){
            shuffle(cycle, cycle+size, g);
            int next;
            cost = 0;
            for (int i = 0; i<size; i++){
                next = (i+1)%size;
                cost += adjacencyMatrix[cycle[i]][cycle[next]];
            }
            if(cost<min){
                min = cost;
            }
        }

        delete[] cycle;

        return min;
    }

    int MST(int start){
        int edgeCount=0;
        nodes[start].visited=true;
        int wage = 0;
        while(edgeCount<size-1){
            //printf("MST:%f\n",float(edgeCount)/float(size-1));
            int min = INT32_MAX;
            int a=-1;
            int b=-1;
            for (int i=0; i < size; i++){
                for(int j=0; j<size; j++)
                    if(adjacencyMatrix[i][j]<min&&isValidEdge(i,j)){
                        min = adjacencyMatrix[i][j];
                        a = i;
                        b = j;
                    }
            }
            if(a!=-1&&b!=-1){
                wage +=min;
                edgeCount++;
                //printf("%d, ", edgeCount);
                nodes[a].visited=true;
                nodes[b].visited=true;
                tree->addNode(a,b);
            }
        }
        //printf("MST:%f\nINITIALIZED",float(edgeCount)/float(size-1));
        return wage;
    }

    int wholeCost(vector<int> cycle){
        int cost = 0;
        for(int i=0; i<size; i++){
            int next = (i+1)%size;
            cost += this->adjacencyMatrix[cycle[i]][cycle[next]];
            //printf("%d ", cycle[i]);
        }
        return cost;
    }

    void swapElementsInVector(std::vector<int>& vec, int index1, int index2) {
        int temporary = vec[index1];
        vec[index1] = vec[index2];
        vec[index2] = temporary;
    }

    int partialCost(int a, int b,vector<int>& cycle){
        int la, lb, ra, rb;
        if (a==0){
            la = size-1;
        }else{
            la = a-1;
        }
        if (b==0){
            lb=size-1;
        }else{
            lb=b-1;
        }
        ra = (a+1)%size;
        rb = (b+1)%size;
        int cost = adjacencyMatrix[cycle[la]][cycle[a]]+ adjacencyMatrix[cycle[a]][cycle[ra]]\
            + adjacencyMatrix[cycle[lb]][cycle[b]]+adjacencyMatrix[cycle[b]][cycle[rb]];
        return cost;
    }

    void swapInts(int* a, int* b){
        int temp = *a;
        *a = *b;
        *b=temp;
    }
    
    double boltzman(int currentEnergy, int neighbourEnergy, double temp){
        return exp((currentEnergy-neighbourEnergy)/temp);
    }
    vector<int> annealing(vector<int> cycle, double temp, double deltaT, int eraIters, bool printCurrent){
        vector<int> current = cycle;
        int bestCost = wholeCost(cycle);
        vector<int> bestSol = cycle;
        uniform_real_distribution<double> annealingDistro;
        uniform_int_distribution<int> distribution(0, cycle.size() - 1);
        int iterations = 0;
        while(temp>1){
            iterations++;
            int index1=distribution(g);
            int index2;
            do{
                index2=distribution(g);
            }while(index1==index2);
            if(index1>index2){
                swapInts(&index1, &index2);
            }

            int partCostBefore = partialCost(index1, index2, current);
            swapElementsInVector(current, index1, index2);
            int partCostAfter = partialCost(index1, index2, current);
            swapElementsInVector(current, index1, index2);
            if(partCostAfter<partCostBefore){
                reverse(current.begin() + index1, current.begin() + index2+1);

                if(wholeCost(current)<wholeCost(bestSol)){
                    //printf("ZOBA");
                    bestSol=current;
                }
                
                // reverse(bestSol.begin() + index1, bestSol.begin() + index2+1);
            }else{
                if(boltzman(partCostBefore, partCostAfter, temp)>annealingDistro(g)){
                    reverse(current.begin() + index1, current.begin() + index2+1);
                }
            }
            if(printCurrent==true && iterations%10000==00){
                printf("%d\n", wholeCost(current));
            }
            if(iterations%eraIters==0){
                temp*=deltaT;
            }
        }
        if(!printCurrent){
            //printf("bestSol = %d, current = %d\n", wholeCost(bestSol), wholeCost(current));
        }
        return bestSol;
    }

    vector<int> localSearch(vector<int> cycle, int &steps){
        int cost = this->wholeCost(cycle);
        vector<int> tempSolution = cycle;
        steps=0;
        while(true){
            int smallestDeltaCost = 0;
            int bestA, bestB;
            for (int i=0; i<size-3; i++){
                //printf("i %d\n", i);
                for(int j=i+3;j<size;j++){
                    //printf("j %d\n", j);
                    int jN=(j+1)%size;
                    int iN = i+1;
                    int jP = j-1;
                    int iP=i-1;
                    int cost;
                    int newCost;
                    if(i==0){
                        iP = size-1;
                    }
                    cost = adjacencyMatrix[tempSolution[iP]][tempSolution[i]]+adjacencyMatrix[tempSolution[i]][tempSolution[iN]]+\
                        adjacencyMatrix[tempSolution[jP]][tempSolution[j]]+adjacencyMatrix[tempSolution[j]][tempSolution[jN]];
                    
                    int temp = tempSolution[i];
                    tempSolution[i]=tempSolution[j];
                    tempSolution[j]=temp;
                    temp = tempSolution[iN];
                    tempSolution[iN]=tempSolution[jP];
                    tempSolution[jP]=temp;

                    newCost = adjacencyMatrix[tempSolution[iP]][tempSolution[i]]+adjacencyMatrix[tempSolution[i]][tempSolution[iN]]+\
                        adjacencyMatrix[tempSolution[jP]][tempSolution[j]]+adjacencyMatrix[tempSolution[j]][tempSolution[jN]];

                    temp = tempSolution[i];
                    tempSolution[i]=tempSolution[j];
                    tempSolution[j]=temp;
                    temp = tempSolution[iN];
                    tempSolution[iN]=tempSolution[jP];
                    tempSolution[jP]=temp;

                    int deltaCost = newCost-cost;
                    if(smallestDeltaCost>deltaCost){
                        smallestDeltaCost = deltaCost;
                        bestA = i;
                        bestB = j;
                    }
                }
            }
            //printf("%d\n", steps);
            if(smallestDeltaCost<0){
                // printf("%d, %d, %d\n", smallestDeltaCost, bestA, bestB);
                reverse(tempSolution.begin() + bestA, tempSolution.begin() + bestB+1);
                // for(int i=0; i<size;i++){
                //     printf("%d ", tempSolution[i]);
                // }
                // printf("\n");
            }else{
                break;
            }
            if(steps%100==0){
                //printf("step mod 100 = %d\n", steps);
            }
            steps++;
        }
        //printf("curren sol: %d\n", wholeCost(tempSolution));
        return tempSolution;
    }

    ////////////////////////////\
    random taboo\
    /////////////////////////
    vector<int> tabuRandom(vector<int> cycle, int tabuSize, int stagnationCount){
        indexPair taboo[tabuSize];
        for(int i=0; i<tabuSize;i++){
            indexPair pair;
            pair.first=-1;
            pair.second=-1;
            taboo[i]=pair;
        }
        vector<int> current = cycle;
        int bestCost = wholeCost(cycle);
        vector<int> bestSol = cycle;
        int iterations = 0;
        int taboCounter = 0;
        int noUpdateCounter=0;
        uniform_int_distribution<int> distribution(0, cycle.size() - 1);
        while (true){
            int bestDelta = INT32_MAX;
            int bestI=-1;
            int bestJ=-1;
            for(int i=0; i<this->size; i++){
                int index1=distribution(g);
                int index2;
                do{
                    index2=distribution(g);
                }while(index1==index2);
                if(index1>index2){
                    swapInts(&index1, &index2);
                }
                int partCostBefore = partialCost(index1, index2, current);
                swapElementsInVector(current, index1, index2);
                int partCostAfter = partialCost(index1, index2, current);
                swapElementsInVector(current, index1, index2);
                int delta = partCostAfter-partCostBefore;
                if(bestDelta>delta){
                    bool inTaboo=false;
                    for(int k=0; k<tabuSize; k++){
                        if(taboo[k].first==index1&&taboo[k].second==index2 || taboo[k].first==index2&&taboo[k].second==index1){
                            inTaboo=true;
                            break;
                        }
                    }
                    if(!inTaboo){
                        bestDelta=delta;
                        bestI=index1;
                        bestJ=index2;
                    }
                }
            }
            taboo[taboCounter].first=bestI;
            taboo[taboCounter].second=bestJ;
            taboCounter++;
            taboCounter%=tabuSize;
            reverse(current.begin() + bestI, current.begin() + bestJ+1);
            if(wholeCost(current)<wholeCost(bestSol)){
                bestSol=current;
                noUpdateCounter=0;
            }else{
                noUpdateCounter++;
            }
            if(noUpdateCounter>=stagnationCount){
                break;
            }
        }
        //printf("current: %d, best: %d\n", wholeCost(current), wholeCost(bestSol));
        return bestSol;
    }

};

int main(int argc, char* argv[]){
    int startMST=0; int startDFS=0;
    if(argc<1){
        printf("Nalezy podac zrodlo grafu");
        return -1;
    }
    ifstream file(argv[1]);

    vector<int> ids;
    vector<int> xs;
    vector<int> ys;
    string line;
    int size = 0;
    std::string input = argv[1];
    size_t dot = input.find_last_of('.');
    size_t slash = input.find_last_of('/');
    std::string toPath = input.substr(slash, dot);

    if (file.is_open()) {
        //printf("Graf otwarty\n");
        while (getline(file, line)) {
            istringstream stream(line);
            int id, x, y;

            if (stream>>id>>x>>y){
                size++;
                ids.push_back(id);
                xs.push_back(x);
                ys.push_back(y);
            }
        }
    }else{
        printf("Nalezy podac prawidlowe zrodlo grafu przy wywolaniu");
        return -1;
    }


    Graph *graph = new Graph(ids.size());

    for(int i=0; i<ids.size(); i++){
        graph->loadNode(xs[i],ys[i],ids[i]-1);
    }
    //printf("Inicjalizacja grafu\n");
    graph->prepareGraph();
    //printf("Graf zinicjalizowany\nObliczanie MST\n");
    //cout<<toPath;

    vector<int> randCycle;
    for(int i =0; i<graph->size; i++){
        randCycle.push_back(i);
    }
    shuffle(begin(randCycle), end(randCycle), g);


    vector<int> bestCycleAnnealin;
    int bestCostAnnealing = INT32_MAX;

    vector<int> bestCycleTabu;
    int bestCostTabu = INT32_MAX;
    unsigned long long avgCostTabu = 0;
    unsigned long long avgCostAnnealing = 0;

    auto total_durationTabu = std::chrono::milliseconds(0);
    auto total_durationAnnealing = std::chrono::milliseconds(0);

    vector<int> cycle=randCycle;
    int iterations = 100;

    for(int i=0; i<iterations;i++){
        shuffle(begin(randCycle), end(randCycle), g);

        auto start_time = std::chrono::high_resolution_clock::now();
        cycle = graph->tabuRandom(randCycle, 7, 1000);
        int costTabu = graph->wholeCost(cycle);
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

        avgCostTabu+=costTabu;
        total_durationTabu+=duration;

        if(costTabu<bestCostTabu){
            bestCostTabu = costTabu;
            bestCycleTabu = cycle;
        }

        start_time = std::chrono::high_resolution_clock::now();
        cycle = graph->annealing(randCycle, 1000, 0.9999, 100, false);

        int costAnnealing = graph->wholeCost(cycle);
        end_time = std::chrono::high_resolution_clock::now();

        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        total_durationAnnealing+=duration;

        if(costAnnealing<bestCostAnnealing){
            bestCostAnnealing = costAnnealing;
            bestCycleAnnealin = cycle;
        }

        avgCostAnnealing +=costAnnealing;
    }
    int steps = 0;
    cycle = graph->localSearch(randCycle, steps);

    int costLocal =graph->wholeCost(cycle);
    avgCostAnnealing/=iterations;
    avgCostTabu/=iterations;
    total_durationAnnealing/=iterations;
    total_durationTabu/=iterations;
    
    cout<<total_durationAnnealing.count()<<endl;
    cout<<total_durationTabu.count()<<endl;
    cout<<toPath<<"&"<<bestCostTabu<<"&"<<avgCostTabu<<"&"<<bestCostAnnealing<<"&"<<avgCostAnnealing<<"&"<<costLocal<<"\\\\ \\hline"<<endl;
    return 0;
}


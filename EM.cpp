#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <cstdlib>

// soft update probability :: add new variable in unknowns (weights)

// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory
using namespace std;

// Our graph consists of a list of nodes where each node is represented as follows:
class Graph_Node{

private:
	string Node_Name;  // Variable name 
	vector<int> Children; // Children of a particular node - these are index of nodes in graph.
	vector<string> Parents; // Parents of a particular node- note these are names of parents
	vector<int> Parents_vals;
	int nvalues;  // Number of categories a variable represented by this node can take
	vector<float> CPT; // conditional probability table as a 1-d array . Look for BIF format to understand its meaning
public:
	vector<int> count;
	vector<int> offset;
	vector<string> values; // Categories of possible values

	// Constructor- a node is initialised with its name and its categories
    Graph_Node(string name,int n,vector<string> vals)
	{
		Node_Name=name;
	
		nvalues=n;
		values=vals;
		

	}
	void inc_count(int pos){
		count[pos]++;
	}
	void dec_count(int pos){
		count[pos] = count[pos]-1;
	}
	void set_offset(vector<int> new_offset){
		offset.clear();
		offset = new_offset;
	}
	vector<int> get_offset(){
		return offset;
	}
	int get_offset(int n){
		return offset[n];
	}
	string get_name()
	{
		return Node_Name;
	}
	vector<int> get_children()
	{
		return Children;
	}
	vector<int> get_Parents()
	{
		return Parents_vals;
	}
	vector<float> get_CPT()
	{
		return CPT;
	}
	float get_CPT(int pos){
		return CPT[pos];
	}
	vector<int> get_count(){
		return count;
	}
	int get_count(int pos){
		return count[pos];
	}
	void count_init(int n){
		count.assign(n, 1);
		return;
	}
	int get_nvalues()
	{
		return nvalues;
	}
	vector<string> get_values()
	{
		return values;
	}
	void set_CPT(vector<float> new_CPT)
	{
		CPT.clear();
		CPT=new_CPT;
	}
    void set_Parents(vector<string> Parent_Nodes)
    {
        Parents.clear();
        Parents=Parent_Nodes;
    }
    void set_Parents_Vals(vector<int> vals){
    	Parents_vals.clear();
    	Parents_vals = vals;
    }
    // add another node in a graph as a child of this node
    int add_child(int new_child_index )
    {
        for(int i=0;i<Children.size();i++)
        {
            if(Children[i]==new_child_index)
                return 0;
        }
        Children.push_back(new_child_index);
        return 1;
    }



};


 // The whole network represted as a list of nodes
class network{
public:

	vector<Graph_Node> Pres_Graph;

	int addNode(Graph_Node node)
	{
		Pres_Graph.push_back(node);
		return 0;
	}
    
    
	int netSize()
	{
		return Pres_Graph.size();
	}
    // get the index of node with a given name
    int get_index(string val_name)
    {
        vector<Graph_Node>::iterator listIt;
        int count=0;
        for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
        {
            if(listIt->get_name().compare(val_name)==0)
                return count;
            count++;
        }
        return -1;
    }
// get the node at nth index
    Graph_Node* get_nth_node(int n)
    {
    	if ( n<Pres_Graph.size() ){
    		return &Pres_Graph[n];
    	}
    	return NULL;
        // vector<Graph_Node>::iterator listIt;
        // int count=0;
        // for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
        // {
        //     if(count==n)
        //         return listIt;
        //     count++;
        // }
        // return listIt; 
    }
    //get the iterator of a node with a given name
    vector<Graph_Node>::iterator search_node(string val_name)
    {
        vector<Graph_Node>::iterator listIt;
        for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
        {
            if(listIt->get_name().compare(val_name)==0)
                return listIt;
        }
    	cout<<"node not found\n";
        return listIt;
    }
	

};

network read_network(string fileName)
{
	network Alarm;
	string line;
	int find=0;
  	ifstream myfile(fileName); 
  	string temp;
  	string name;
  	vector<string> values;
  	vector<int> vals;
    if (myfile.is_open())
    {
    	while (! myfile.eof() )
    	{
    		stringstream ss;
      		getline (myfile,line);
      		
      		
      		ss.str(line);
     		ss>>temp;
     		
     		
     		if(temp.compare("variable")==0)
     		{
                    
     				ss>>name;
     				getline (myfile,line);
                   
     				stringstream ss2;
     				ss2.str(line);
     				for(int i=0;i<4;i++)
     				{
     					
     					ss2>>temp;
     					
     					
     				}
     				values.clear();
     				while(temp.compare("};")!=0)
     				{
     					values.push_back(temp);
     					
     					ss2>>temp;
    				}
     				Graph_Node new_node(name,values.size(),values);
     				int pos=Alarm.addNode(new_node);

     				
     		}
     		else if(temp.compare("probability")==0)
     		{
                    
     				ss>>temp;
     				ss>>temp;
     				
                    vector<Graph_Node>::iterator listIt;
                    vector<Graph_Node>::iterator listIt1;
     				listIt=Alarm.search_node(temp);
                    int index=Alarm.get_index(temp);
                    ss>>temp;
                    values.clear();
                    vals.clear();
     				while(temp.compare(")")!=0)
     				{
                        listIt1=Alarm.search_node(temp);
                        listIt1->add_child(index);
     					values.push_back(temp);
     					vals.push_back(Alarm.get_index(temp));
     					ss>>temp;

    				}
                    listIt->set_Parents(values);
                    listIt->set_Parents_Vals(vals);
    				getline (myfile,line);
     				stringstream ss2;
                    
     				ss2.str(line);
     				ss2>> temp;
                    
     				ss2>> temp;
                    
     				vector<float> curr_CPT;
                    string::size_type sz;
     				while(temp.compare(";")!=0)
     				{
                        
     					curr_CPT.push_back(atof(temp.c_str()));
     					
     					ss2>>temp;
                     
    				}
                    
                    listIt->set_CPT(curr_CPT);
                    listIt->count_init(curr_CPT.size());
     		}
            else
            {
                
            }
     	
    	}

  	}
  	
  	myfile.close();
  	return Alarm;
}

void initialize_network(network &Alarm){
	vector<Graph_Node>::iterator listIt;
	int j=0;
    for(listIt=Alarm.Pres_Graph.begin(); listIt!=Alarm.Pres_Graph.end(); listIt++)
    {
        vector<float> new_CPT = listIt->get_CPT();
        for(int i=0; i<new_CPT.size(); i++){
        	new_CPT[i] = (float)rand()/RAND_MAX;
        }
        listIt->set_CPT(new_CPT);

        vector<int> parents = listIt->get_Parents();
        int pSize[parents.size()];
        int mul=1;
        for(int i=0; i<parents.size(); i++){
        	pSize[i] = Alarm.get_nth_node(parents[i])->get_nvalues();
        	mul *= pSize[i];
        }
        vector<int> offset;
        for(int i=0 ; i<parents.size()+1; i++){
        	offset.push_back(mul);
        	if( i< parents.size())	mul /= pSize[i];
        }
        listIt->set_offset(offset);
    }
}

void Write_file(network &Alarm, string fileName){
	ofstream outFile("solved_alarm.bif");

	string line;
	ifstream myfile(fileName); 
  	string temp;
  	string name;
  	vector<string> values;
  	
    if (myfile.is_open())
    {
    	while ( !myfile.eof() )
    	{
    		stringstream ss;
      		getline (myfile,line);
      		
      		
      		ss.str(line);
     		ss>>temp;
     		if (line.compare("")==0){
     			continue;
     		}
     		else if(temp.compare("probability")==0)
     		{
     				line += "\n";
                    outFile << line;
     				ss>>temp;
     				ss>>temp;
     				
                    vector<Graph_Node>::iterator listIt;
                    listIt=Alarm.search_node(temp);
                    vector<float> CPT = listIt->get_CPT();
                    line = "\ttable";
    				for(int i=0; i<CPT.size(); i++){
    					line += " ";
    					line += to_string(CPT[i]); 
    				}
    				line += " ;\n";
    				outFile << line;    				
     		}
            else if ((temp.compare("table") != 0) && (temp.compare("probability") != 0))
            {
            	line += "\n";
                outFile << line;
            }
    	}

  		myfile.close();
  	}
  	outFile.close();
}

// node val curr_offset child_offset chil
vector<vector<int>> unknowns;

int assign_val(Graph_Node node, vector<int> parents_val){
	int offset = 0;
	for(int i=0; i<parents_val.size(); i++){
		offset += node.offset[i+1] * parents_val[i];
	}
	vector<float> CPT = node.get_CPT();
	float m = CPT[offset];
	int index = 0;
	for(int i=1; i<node.get_nvalues(); i++){
		if( CPT[i*node.offset[0] + offset ] > m ){
			m = CPT[i*node.offset[0] + offset ];
			index = i;
		}
	}
	return index;
}

int get_Updateoffset(network &Alarm, int n, int p, int *assigned_val, int &p_pos){
	Graph_Node *node = Alarm.get_nth_node(n);
	vector<int> parents = node->get_Parents();

	vector<int> parents_val;
	for(int i=0; i<parents.size(); i++){
		parents_val.push_back(assigned_val[parents[i]]);
	}
	int offset = assigned_val[n]*(node->get_offset(0));
	for(int i=0; i < parents_val.size(); i++){
		if(parents[i] != p){
			offset += parents_val[i]*(node->get_offset(i+1));
		}
		else{
			p_pos = node->get_offset(i+1);
		}
	}
	return offset;
}

void initialize_unknowns(string fileName, network &Alarm){
	ifstream dFile(fileName);
	if (!dFile.is_open()){
		cout << "could not open file" <<endl;
		return;
	}
	int n = Alarm.netSize();
	int assigned_val[n];
	while ( !dFile.eof() )
	{
		string line = "";
		getline(dFile, line);
		stringstream ss(line);
		string temp;
		int store = -1;
		for(int i=0; i<n; i++){
			ss >> temp;
			if (temp == "\"?\""){
				store = i;
			}
			else{
				Graph_Node *temp_graph = Alarm.get_nth_node(i);
				for(int j=0; j<temp_graph->values.size(); j++){
					if(temp_graph->values[j].compare(temp) == 0){
						assigned_val[i] = j;
						break;
					}
				}
			}
		}
		// only if store >= 0
		if(store >= 0){
				Graph_Node node = Alarm.Pres_Graph[store];
				vector<int> parents = node.get_Parents();

				vector<int> parents_val;
				for(int i=0; i<parents.size(); i++){
					parents_val.push_back(assigned_val[parents[i]]);
				}
				assigned_val[store] = assign_val(node, parents_val);

				vector<int> unknown_val;
				unknown_val.push_back(store);
				unknown_val.push_back(assigned_val[store]);
				vector<int> offset = node.get_offset();
				int val=0;
				for(int i=0; i<parents_val.size(); i++){
					val += offset[i+1]*parents_val[i];
				}
				unknown_val.push_back(val);
				vector<int> children = node.get_children();
				for(int i=0; i<children.size(); i++){
					int p_pos;
					unknown_val.push_back(get_Updateoffset(Alarm, children[i], store, assigned_val, p_pos));
					unknown_val.push_back(p_pos);
				}
				unknowns.push_back(unknown_val);
			}

		// increase node count
		for(int i=0; i<n; i++){
			int unused = 0;
			int pos = get_Updateoffset(Alarm, i, -1, assigned_val, unused);
			Alarm.Pres_Graph[i].inc_count(pos);
		}
	}
	dFile.close();
}

int EM(network &Alarm){
	int update_cnt = 0;
	// calculate probability
	vector<Graph_Node>::iterator listIt;
	for(listIt = Alarm.Pres_Graph.begin(); listIt!= Alarm.Pres_Graph.end(); listIt++){
		vector<int> cnt = listIt->get_count();
		vector<float> CPT = listIt->get_CPT();
		int offset = listIt->get_offset(0);
		for(int j=0; j<offset; j++){
			int i = 0;
			int s = 0;
			while(i*offset < cnt.size()){
				s += cnt[j + i*offset];
				i++;
			}
			i = 0;
			while(i*offset < cnt.size()){
				CPT[j + i*offset] = (float)cnt[j+i*offset]/s;
				i++;
			}
		}
		listIt->set_CPT(CPT);
	}

	// update counts
	for(int i=0; i<unknowns.size(); i++){
		int node = unknowns[i][0];
		vector<int> children = Alarm.Pres_Graph[node].get_children();
		int offset = Alarm.Pres_Graph[node].get_offset(0);

		int old_val = unknowns[i][1];
		// assign new value to unknown
		int nvalues = Alarm.Pres_Graph[node].get_nvalues();
		int new_val = 0; // value at which probability is maximum
		float prev_prob = 0;
		for(int k=0; k < nvalues; k++){
			float prob = Alarm.Pres_Graph[node].get_CPT(k*offset + unknowns[i][2]);
			for(int j=0; j<children.size(); j++){
				prob *= Alarm.Pres_Graph[children[j]].get_CPT(unknowns[i][3+j*2] + k*unknowns[i][3+j*2+1]);
			}
			if(prob > prev_prob){
				prev_prob = prob;
				new_val = k;
			}
		}

		// decrease count
		Alarm.Pres_Graph[node].dec_count(old_val*offset + unknowns[i][2]);
		// if(Alarm.Pres_Graph[node].get_count(old_val*offset + unknowns[i][2]) < 0) { cout << "less than 0: "<< node << " " << old_val << " " << Alarm.Pres_Graph[node].get_count(old_val*offset + unknowns[i][2]) << " : unknown " << i <<endl;}
		for(int j=0; j<children.size(); j++){
			Alarm.Pres_Graph[children[j]].dec_count(unknowns[i][3+j*2] + old_val*unknowns[i][3+j*2+1]);
			// if(Alarm.Pres_Graph[children[j]].get_count(unknowns[i][3+j*2] + old_val*unknowns[i][3+j*2+1]) < 0){ cout << "less than 0: "<< children[i] <<endl;}
		}

		// increase new count
		Alarm.Pres_Graph[node].inc_count(new_val*offset + unknowns[i][2]);
		for(int j=0; j<children.size(); j++){
			Alarm.Pres_Graph[children[j]].inc_count(unknowns[i][3+j*2] + new_val*unknowns[i][3+j*2+1]);
		}

		// update unknowns
		unknowns[i][1] = new_val;
		if (new_val != old_val) update_cnt++;
	}
	return update_cnt;
}

void update_prob(network &Alarm){
	vector<Graph_Node>::iterator listIt;
	for(listIt = Alarm.Pres_Graph.begin(); listIt!= Alarm.Pres_Graph.end(); listIt++){
		vector<int> cnt = listIt->get_count();
		vector<float> CPT = listIt->get_CPT();
		int offset = listIt->get_offset(0);
		for(int j=0; j<offset; j++){
			int i = 0;
			int s = 0;
			while(i*offset < cnt.size()){
				s += cnt[j + i*offset];
				i++;
			}
			i = 0;
			while(i*offset < cnt.size()){
				CPT[j + i*offset] = (float)cnt[j+i*offset]/s;
				i++;
			}
		}
		listIt->set_CPT(CPT);
	}
	return;
}

int main(int argc, char const *argv[])
{
	network Alarm;
	network temp_Alarm;
	string dataFile(argv[2]);
	string inpFile(argv[1]);
	Alarm=read_network(inpFile);
    initialize_network(Alarm);
    initialize_unknowns(dataFile, Alarm);
    temp_Alarm = Alarm;

    int update_cnt;
    for (int i=0; i<40; i++){
    	// cout << i << ": ";
        update_cnt = EM(temp_Alarm);
        // cout << update_cnt << endl;
        if(update_cnt == 0) {
        	Alarm = temp_Alarm;
        	initialize_unknowns(dataFile, temp_Alarm);
        	// break;
        }
    }
    update_prob(Alarm);
    Write_file(Alarm, inpFile);
    return 0;
}
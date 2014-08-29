#include <iostream>
#include <cmath>
#include <exception>
using namespace std;

struct polarNum;
class connection;
class node;

polarNum polarParse(float,float);

struct polarNum 
{
	float l;
	float a;
	polarNum operator+(polarNum &b)
	{
		return polarParse(this->l*cos(this->a)+b.l*cos(b.a),this->l*sin(this->a)+b.l*sin(b.a));
	}
	polarNum operator-(polarNum &b)
	{
		return polarParse(this->l*cos(this->a)-b.l*cos(b.a),this->l*sin(this->a)-b.l*sin(b.a));
	}
	polarNum operator/(polarNum &b)
	{
		polarNum temp;
		temp.l = this->l/b.l;
		temp.a = this->a-b.a;
		return temp;
	}
	polarNum operator*(polarNum &b)
	{
		polarNum temp;
		temp.l = this->l*b.l;
		temp.a = this->a+b.a;
		return temp;
	}
};

class node
{
private:
	int connectionCount;
	int connectionCounter;
	int selfNum;
	connection **connectionArray;

public:
	node();
	~node();
	void setConnectionCount(int,int);
	void addConnectionDetails(connection*);
	connection* getNotTraversedConnection();
	connection* getConnectionTo(node*);
	connection* getNextConnection();
	int getSelfNum();
	void resetCounter();
};

class connection
{
private:
	node *fromNode;
	node *toNode;
	int traversed;
	polarNum impedence;
	polarNum voltage;
	polarNum current;
public:
	connection();
	~connection();
	void setDetails(node*,node* /*,polarNum,polarNum,polarNum*/);
	void printData();
	void setTraversed();
	void setUntraversed();
	bool isTraversed();
	node* getOtherNode(node*);
};

struct currentPath
{
	node *through;
	connection *from;
	connection *to;
	currentPath *prev;
	currentPath *next;
	currentPath *suppose;
	static int pathCountActual;
	static int pathCountSuppose;
	currentPath(node*,currentPath*);
	static currentPath* make(node*,currentPath*);
	static currentPath* start;
	static currentPath* startNewTraversePath(node*);
	currentPath* getNextNode();
	bool hasOccured(node*);
};

node::node()
{
	connectionCounter = 0;
}

void node::setConnectionCount(int n,int a)
{
	this->connectionArray = new connection*[n];
	this->connectionCount = n;
	this->selfNum = a;
}

void node::addConnectionDetails(connection *_new)
{
	connectionArray[connectionCounter] = _new;
	connectionCounter++;
	if(connectionCounter == connectionCount)
		connectionCounter = 0;
}

void node::resetCounter()
{
	connectionCounter = 0;
}

connection* node::getNotTraversedConnection()
{
	for (int i = 0; i < connectionCount; ++i)
	{
		if (!connectionArray[i]->isTraversed())
		{
			return connectionArray[i];
		}
	}
	return NULL;
}

connection* node::getConnectionTo(node* _node)
{
	for (int i = 0; i < connectionCount; ++i)
	{
		if(connectionArray[i]->getOtherNode(this) == _node)
			return connectionArray[i];
	}
	return NULL;
}

connection* node::getNextConnection()
{
	if(connectionCounter==connectionCount)
		return NULL;
	return connectionArray[connectionCounter++];
}

int node::getSelfNum()
{
	return this->selfNum;
}

connection::connection()
{
	traversed = 0;
}

void connection::setDetails(node *prev, node *next/*,polarNum imp = 0,polarNum volt = 0,polarNum cur = 0*/)
{
	this->fromNode = prev;
	this->toNode = next;
	//this->impedence = imp;
	//this->voltage = volt;
	//this->current = cur;
}

void connection::printData()
{
	cout<<impedence.l<<"  "<<impedence.a/M_PI*180<<" : "<<voltage.l<<"  "<<voltage.a/M_PI*180<<" : "<<current.l<<"  "<<current.a/M_PI*180;
}

void connection::setTraversed()
{
	this->traversed++;
}

void connection::setUntraversed()
{
	this->traversed--;
}

bool connection::isTraversed()
{
	return this->traversed==0?false:true;
}

node* connection::getOtherNode(node* _node)
{
	return _node==toNode?fromNode:toNode;
}

currentPath* currentPath::start = NULL;
int currentPath::pathCountActual = 0;
int currentPath::pathCountSuppose = 0;

currentPath::currentPath(node* _through,currentPath* _prev)
{
	this->through = _through;
	this->prev = _prev;
	if(_prev == NULL)
	{
		start = this;
		pathCountActual = 0;
		pathCountSuppose = 0;
	}
	else
	{
		_prev->next = this;
		this->from = _prev->to;
	}
}

currentPath* currentPath::make(node* _through,currentPath* _prev = NULL)
{
	return new currentPath(_through,_prev);
}

bool currentPath::hasOccured(node* _node)
{
	currentPath *temp = this;
	while(temp!=NULL){
		if(temp->through == _node)
			return true;
		temp = temp->prev;
	}
	return false;
}

currentPath* currentPath::getNextNode()
{
	connection* temp;
	currentPath* curTemp;
	//cout<<this->through->getSelfNum()<<" > ";
	if(this==start)
	{
		if((temp=this->through->getNotTraversedConnection())!=NULL)
		{
			this->to = temp;
			temp->setTraversed();
			return currentPath::make(temp->getOtherNode(this->through),this)->getNextNode();
		}
		else
		{
			return NULL;
		}
	}
	if(this->prev->through != start->through &&(temp=this->through->getConnectionTo(start->through))!=NULL)
	{
		this->to = temp;
		temp->setTraversed();
		this->next = start;
		start->prev = this;
		return start;
	}
	this->through->resetCounter();
	while((temp=this->through->getNextConnection())!=NULL)
	{
		if(!hasOccured(temp->getOtherNode(this->through)))
		{
			this->to = temp;
			temp->setTraversed();
			curTemp = currentPath::make(temp->getOtherNode(this->through),this)->getNextNode();
			if(curTemp == NULL){
				this->to = NULL;
				temp->setUntraversed();
				delete this->next;
			}
			else
				return curTemp;
		}
	}
	return NULL;
}


currentPath* currentPath::startNewTraversePath(node* _start)
{
	return currentPath::make(_start)->getNextNode();
}

polarNum polarParse(float a, float b)
{
	polarNum temp;
	temp.l = sqrt(a*a+b*b);
	temp.a = atan(b/a);
	return temp;
}

int main()
{	
	node *nodes;
	connection *connections;
	currentPath *one,*two;
	int nodeCount=0,connectionCount=0,temp;
	
	cout<<"Enter the number of nodes : ";
	cin>>nodeCount;
	nodes = new node[nodeCount];
	
	cout<<"Enter the number of connections per node\n";
	for(int i=0;i<nodeCount;i++)	//setting node deatils
	{
		cout<<"Node #"<<i+1<<" : ";
		cin>>temp;
		nodes[i].setConnectionCount(temp,i);
		connectionCount+=temp;
	}
	connectionCount/=2;
	connections = new connection[connectionCount];

	cout<<"Enter the connection details as 'nodeA nodeB Resistance Inductance Capacitance VoltageSourcePeak VoltageSourcePhase CurrentSourcePeak CurrentSourcePhase'\n";
	for(int i=0;i<connectionCount;i++)	//setting connections details
	{
		int nodeA,nodeB;
		float Resistance,Inductance,Capacitance,VoltageSourcePeak,VoltageSourcePhase,CurrentSourcePeak,CurrentSourcePhase;
		cout<<"Connection #"<<i+1<<" : ";
		cin>>nodeA>>nodeB/*>>Resistance>>Inductance>>Capacitance>>VoltageSourcePeak>>VoltageSourcePhase>>CurrentSourcePeak>>CurrentSourcePhase*/;
		connections[i].setDetails(&nodes[nodeA-1],&nodes[nodeB-1]/*,polarParse(Resistance,Inductance-Capacitance),polarNum({VoltageSourcePeak,VoltageSourcePhase/180*M_PI}),polarNum({CurrentSourcePeak,CurrentSourcePhase/180*M_PI})*/);
		nodes[nodeA-1].addConnectionDetails(&connections[i]);
		nodes[nodeB-1].addConnectionDetails(&connections[i]);
	}

	for (int i = 0; i < nodeCount; ++i)
	{
		while((one = currentPath::startNewTraversePath(&nodes[i]))!=NULL)
		{
			two = one;
			do
			{
				cout<<two->through->getSelfNum()+1<<" > ";
				two = two->next;
			} while (two!=one);
			cout<<two->through->getSelfNum()+1<<endl;
		}
	}

	system("pause");
	system("cls");
	return 0;
}
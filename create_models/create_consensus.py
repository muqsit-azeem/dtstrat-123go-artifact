def get_consensus_str(N):
	s = f"""// COIN FLIPPING PROTOCOL FOR POLYNOMIAL RANDOMIZED CONSENSUS [AH90] 
// gxn/dxp 20/11/00

mdp

// constants
const int N={N};
const int K;
const int range = 2*(K+1)*N;
const int counter_init = (K+1)*N;
const int left = N;
const int right = 2*(K+1)*N - N;

// shared coin
global counter : [0..range] init counter_init;

module process1
	
	// program counter
	pc1 : [0..3];
	// 0 - flip
	// 1 - write 
	// 2 - check
	// 3 - finished
	
	// local coin
	coin1 : [0..1];	

	// flip coin
	[] (pc1=0)  -> 0.5 : (coin1'=0) & (pc1'=1) + 0.5 : (coin1'=1) & (pc1'=1);
	// write tails -1  (reset coin to add regularity)
	[] (pc1=1) & (coin1=0) & (counter>0) -> (counter'=counter-1) & (pc1'=2) & (coin1'=0);
	// write heads +1 (reset coin to add regularity)
	[] (pc1=1) & (coin1=1) & (counter<range) -> (counter'=counter+1) & (pc1'=2) & (coin1'=0);
	// check
	// decide tails
	[] (pc1=2) & (counter<=left) -> (pc1'=3) & (coin1'=0);
	// decide heads
	[] (pc1=2) & (counter>=right) -> (pc1'=3) & (coin1'=1);
	// flip again
	[] (pc1=2) & (counter>left) & (counter<right) -> (pc1'=0);
	// loop (all loop together when done)
	[done] (pc1=3) -> (pc1'=3);

endmodule

// construct remaining processes through renaming"""

	for i in range(1,N):
		s += f"\nmodule process{i+1} = process1[pc1=pc{i+1},coin1=coin{i+1}] endmodule"

	s += """\n
// labels
"""
	s += f"label \"finished\" = "
	for i in range(N):
		s += f"pc{i+1}=3"
		if i != N-1:
			s += " & "
	s += " ;\n"

	s += f"label \"all_coins_equal_0\" = "
	for i in range(N):
		s += f"coin{i+1}=0"
		if i != N-1:
			s += " & "
	s += " ;\n"

	s += f"label \"all_coins_equal_1\" = "
	for i in range(N):
		s += f"coin{i+1}=1"
		if i != N-1:
			s += " & "
	s += " ;\n"

	s += f"label \"agree\" = "
	for i in range(1,N):
		s += f"coin{i}=coin{i+1}"
		if i != N-1:
			s += " & "
	s += " ;\n"

	s+="""
// rewards
rewards "steps"
	true : 1;
endrewards
"""

	# print(s)

	return s

def create_consensus(l):
	for N in l:
		s = get_consensus_str(N)
		f = open(f"consensus.{N}.prism","w")
		f.write(s)
		f.close()

if __name__ == "__main__":
	NList = [2,4,6,8,10] # number of processes (N) you want to create prismfiles for 
	# K will be given as storm arguments
	create_consensus(NList)

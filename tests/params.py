class PARAMS:
    def __init__(self):
        self.N_AFR = 5_000  # african effective population size
        self.N_ANC = 11_000  # ancestral effective population size
        self.N_EUR = 5_000  # european effective population size
        self.N_ND = 1_000  # neanderthal effective population size
        self.N_NND = 10_000  # non-neanderthal effective population size
        self.t_all = 28_000  # time of merging into ANC
        self.t_admix = 2_000  # time of admixture
        self.t_ND = 1_500  # time of neanderthal sample
        self.t_AFR_EUR = 4_000  # time of AFR-EUR merge
        self.TESTCOUNT = 200  # number of tests to run
        self.threshold = 0.8  # threshold for saying sample has admixture
        self.mu = 1.2 * 10 ** -8  # mutation rate
        self.seq_len = 50_000  # sequence length
Params = PARAMS()

# Define parameters

# ------ Gaussian #1 -------

CONSTANT_GAUS1 = 100;
MEAN_GAUS1 = 0;
SIGMA_GAUS1 = 1;

CONSTANT_MAX_GAUS1 = 200;
CONSTANT_MIN_GAUS1 = 0;
MEAN_MAX_GAUS1 = -1;
MEAN_MIN_GAUS1 = 1;
SIGMA_MAX_GAUS1 = 2;
SIGMA_MIN_GAUS1 = 0;

# -------------------------
# ------ Gaussian #2 -------

CONSTANT_GAUS2 = 100;
MEAN_GAUS2 = 0;
SIGMA_GAUS2 = 1;

CONSTANT_MAX_GAUS2 = 200;
CONSTANT_MIN_GAUS2 = 0;
MEAN_MAX_GAUS2 = -1;
MEAN_MIN_GAUS2 = 1;
SIGMA_MAX_GAUS2 = 2;
SIGMA_MIN_GAUS2 = 5;

# -------------------------
# ------ Crystal ball #1 -------

CONSTANT_CRYSTAL1 = 100;
MEAN_CRYSTAL1 = 0;
SIGMA_CRYSTAL = 1;
ALPHA_CRYSTAL = -0.1;
N_CRYSTAL = 3;

CONSTANT_MAX_CRYSTAL1 = 200;
CONSTANT_MIN_CRYSTAL1 = 0;
MEAN_MAX_CRYSTAL1 = -1;
MEAN_MIN_CRYSTAL1 = 1;
SIGMA_MAX_CRYSTAL = 5;
SIGMA_MIN_CRYSTAL = 0;
ALPHA_MAX_CRYSTAL = 0.0;
ALPHA_MIN_CRYSTAL = -0.5;
N_MAX_CRYSTAL = 10;
N_MIN_CRYSTAL = 0;

# ------------------------------
# ------ Crystal ball #2 -------

CONSTANT_CRYSTAL2 = 100;
MEAN_CRYSTAL2 = 0;
SIGMA_CRYSTAL = 1;
ALPHA_CRYSTAL = -0.1;
N_CRYSTAL = 3;

CONSTANT_MAX_CRYSTAL2 = 200;
CONSTANT_MIN_CRYSTAL2 = 0;
MEAN_MAX_CRYSTAL2 = -1;
MEAN_MIN_CRYSTAL2 = 1;
SIGMA_MAX_CRYSTAL = 5;
SIGMA_MIN_CRYSTAL2 = 0;
ALPHA_MAX_CRYSTAL2 = 0.0;
ALPHA_MIN_CRYSTAL2 = -0.5;
N_MAX_CRYSTAL2 = 10;
N_MIN_CRYSTAL2 = 0;

# -----------------------------


#----------------------------------------
# Import and set parameters from the file
#----------------------------------------

def set_params(filename, fitType):
    parFile = open(filename)
    table = csv.reader(parFile, delimiter='\t')
    if (len(table)>2):
       PRINTER.print("Invalid parameter set format!")
       parser.print_help()
       sys.exit()

   table = csv.reader()
   # validate the table

   # set parameters
   row = table[0]

   if (fitType==1 or fitType==3 or fitType==5):
      CONSTANT_GAUS1 = row[0]
      CONSTANT_MIN_GAUS1 = row[1]
      CONSTANT_MAX_GAUS1 = row[2]
      MEAN_GAUS1 = row[3]
      MEAN_MIN_GAUS1 = row[4]
      MEAN_MAX_GAUS1 = row[5]
      SIGMA_GAUS1 = row[6]
      SIGMA_MIN_GAUS1 = row[7]
      SIGMA_MAX_GAUS1 = row[8]

   
   if (fitType==2 or fitType==4):
      CONSTANT_CRYSTAL1 = row[0]
      CONSTANT_MIN_CRYSTAL1 = row[1]
      CONSTANT_MAX_CRYSTAL1 = row[2]
      MEAN_CRYSTAL1 = row[3]
      MEAN_MIN_CRYSTAL1 = row[4]
      MEAN_MAX_CRYSTAL1 = row[5]
      SIGMA_CRYSTAL1 = row[6]
      SIGMA_MIN_CRYSTAL1 = row[7]
      SIGMA_MAX_CRYSTAL1 = row[8]
      ALPHA_CRYSTAL1 = row[9]
      ALPHA_MIN_CRYSTAL1 = row[10]
      ALPHA_MAX_CRYSTAL1 = row[11]
      N_CRYSTAL1 = row[12]
      N_MIN_CRYSTAL1 = row[13]
      N_MAX_CRYSTAL1 = row[14]

   row = table[1]   
   
   if (fitType==5):
      CONSTANT_CRYSTAL2 = row[0]
      CONSTANT_MIN_CRYSTAL2 = row[1]
      CONSTANT_MAX_CRYSTAL2 = row[2]
      MEAN_CRYSTAL2 = row[3]
      MEAN_MIN_CRYSTAL2 = row[4]
      MEAN_MAX_CRYSTAL2 = row[5]
      SIGMA_CRYSTAL2 = row[6]
      SIGMA_MIN_CRYSTAL2 = row[7]
      SIGMA_MAX_CRYSTAL2 = row[8]
      ALPHA_CRYSTAL2 = row[9]
      ALPHA_MIN_CRYSTAL2 = row[10]
      ALPHA_MAX_CRYSTAL2 = row[11]
      N_CRYSTAL2 = row[12]
      N_MIN_CRYSTAL2 = row[13]
      N_MAX_CRYSTAL2 = row[14]

   if (fitType==3):
      CONSTANT_GAUS2 = row[0]
      CONSTANT_MIN_GAUS2 = row[1]
      CONSTANT_MAX_GAUS2 = row[2]
      MEAN_GAUS2 = row[3]
      MEAN_MIN_GAUS2 = row[4]
      MEAN_MAX_GAUS2 = row[5]
      SIGMA_GAUS2 = row[6]
      SIGMA_MIN_GAUS2 = row[7]
      SIGMA_MAX_GAUS2 = row[8]
   
   if (fitType==4):
      CONSTANT_CRYSTAL2 = row[0]
      CONSTANT_MIN_CRYSTAL2 = row[1]
      CONSTANT_MAX_CRYSTAL2 = row[2]
      MEAN_CRYSTAL2 = row[3]
      MEAN_MIN_CRYSTAL2 = row[4]
      MEAN_MAX_CRYSTAL2 = row[5]
      SIGMA_CRYSTAL2 = row[6]
      SIGMA_MIN_CRYSTAL2 = row[1]
      SIGMA_MAX_CRYSTAL2 = row[2]
      ALPHA_CRYSTAL2 = row[9]
      ALPHA_MIN_CRYSTAL2 = row[10]
      ALPHA_MAX_CRYSTAL2 = row[11]
      N_CRYSTAL2 = row[12]
      N_MIN_CRYSTAL2 = row[13]
      N_MAX_CRYSTAL2 = row[14]

#----------------------------------------

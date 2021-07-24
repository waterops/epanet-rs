/*
 ******************************************************************************
 Project:      OWA EPANET
 Version:      2.2
 Module:       types.h
 Description:  symbolic constants and data types used throughout EPANET
 Authors:      see AUTHORS
 Copyright:    see AUTHORS
 License:      see LICENSE
 Last Updated: 07/11/2020
 ******************************************************************************
*/
#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(dead_code)]
use crate::linked_list::List;
use std::collections::HashMap;
use std::fs::File;

/*
----------------------------------------------
   Various constants
----------------------------------------------
*/
const CODEVERSION: i32 = 20200;
const MAGICNUMBER: i32 = 516114521;
const ENGINE_VERSION: i32 = 201;   // Used for binary hydraulics file
const EOFMARK: i32 = 0x1A;  // Use 0x04 for UNIX systems
const MAXTITLE: i32 = 3;        // Max. # title lines
const TITLELEN: i32 = 79;       // Max. # characters in a title line
const MAXID: i32 = 31;       // Max. # characters in ID name
const MAXMSG: i32 = 255;      // Max. # characters in message text
const MAXLINE: i32 = 1024;     // Max. # characters read from input line
const MAXFNAME: i32 = 259;      // Max. # characters in file name
const MAXTOKS: i32 = 40;       // Max. items per line of input
const TRUE: i32 = 1;
const FALSE: i32 = 0;
const FULL: i32 = 2;

const BIG: f64 = 1.0e10;
const TINY: f64 = 1.0e-6;
const MISSING: f64 = -1.0e10;     // Missing value indicator
const DIFFUS: f64 = 1.3e-8;     // Diffusivity of chlorine
                                // @ 20 deg C (sq ft/sec)
const VISCOS: f64 = 1.1e-5;     // Kinematic viscosity of water
                                // @ 20 deg C (sq ft/sec)
const MINPDIFF: f64 =  0.1;        // PDA min. pressure difference (psi or m)
const SEPSTR: &str = " \t\n\r";  // Token separator characters


/*
----------------------------------------------
   Enumerated Data Types
----------------------------------------------
*/
enum ObjectType {
    NODE, 
    LINK,
    TIMEPAT, 
    CURVE, 
    CONTROL, 
    RULE
}

enum NodeType {
    JUNCTION, 
    RESERVOIR, 
    TANK
}

enum LinkType {
    CVPIPE,        // pipe with check valve
    PIPE,          // pipe
    PUMP,          // pump
    PRV,           // pressure reducing valve
    PSV,           // pressure sustaining valve
    PBV,           // pressure breaker valve
    FCV,           // flow control valve
    TCV,           // throttle control valve
    GPV            // general purpose valve
}

enum HydFileType {
    USE,           // use hydraulics file from previous run
    SAVE,          // save hydraulics file after current run
    SCRATCH        // use temporary hydraulics file
}

enum QualType {
    NONE,          // no quality analysis
    CHEM,          // analyze a chemical
    AGE,           // analyze water age
    TRACE          // trace % of flow from a source
}

enum CurveType {
    VOLUME_CURVE,  // volume curve
    PUMP_CURVE,    // pump curve
    EFFIC_CURVE,   // efficiency curve
    HLOSS_CURVE,   // head loss curve
    GENERIC_CURVE  // generic curve
}

enum PumpType {
    CONST_HP,      // constant horsepower
    POWER_FUNC,    // power function
    CUSTOM,        // user-defined custom curve
    NOCURVE
}

enum SourceType {
    CONCEN,        // inflow concentration
    MASS,          // mass inflow booster
    SETPOINT,      // setpoint booster
    FLOWPACED      // flow paced booster
}

enum ControlType {
    LOWLEVEL,      // act when grade below set level
    HILEVEL,       // act when grade above set level
    TIMER,         // act when set time reached
    TIMEOFDAY      // act when time of day occurs
}

enum StatusType {
    XHEAD,         // pump cannot deliver head (closed)
    TEMPCLOSED,    // temporarily closed
    CLOSED,        // closed
    OPEN,          // open
    ACTIVE,        // valve active (partially open
    XFLOW,         // pump exceeds maximum flow
    XFCV,          // FCV cannot supply flow
    XPRESSURE,     // valve cannot supply pressure
    FILLING,       // tank filling
    EMPTYING,      // tank emptying
    OVERFLOWING    // tank overflowing
}

enum HeadLossType {
    HW,            // Hazen-Williams
    DW,            // Darcy-Weisbach
    CM             // Chezy-Manning
}

enum UnitsType {
    US,            // US
    SI             // SI (metric) 
}

enum FlowUnitsType {
    CFS,           // cubic feet per second
    GPM,           // gallons per minute
    MGD,           // million gallons per day
    IMGD,          // imperial million gal. per day
    AFD,           // acre-feet per day
    LPS,           // liters per second
    LPM,           // liters per minute
    MLD,           // megaliters per day
    CMH,           // cubic meters per hour
    CMD            // cubic meters per day  
}

enum PressureUnitsType {
    PSI,           // pounds per square inch
    KPA,           // kiloPascals
    METERS         // meters
}

enum RangeType {
    LOW,           // lower limit
    HI,            // upper limit
    PREC           // precision
}

enum MixType {
    MIX1,          // complete mix model
    MIX2,          // 2-compartment model
    FIFO,          // first in, first out model
    LIFO           // last in, first out model
}

enum StatisticType {
    SERIES,        // point time series
    AVG,           // time-averages
    MIN,           // minimum values
    MAX,           // maximum values
    RANGE          // max - min values
}

enum FieldType {
    ELEV = 0,      // nodal elevation
    DEMAND,        // nodal demand flow
    HEAD,          // nodal hydraulic head
    PRESSURE,      // nodal pressure
    QUALITY,       // nodal water quality
  
    LENGTH,        // link length
    DIAM,          // link diameter
    FLOW,          // link flow rate
    VELOCITY,      // link flow velocity
    HEADLOSS,      // link head loss
    LINKQUAL,      // avg. water quality in link
    STATUS,        // link status
    SETTING,       // pump/valve setting
    REACTRATE,     // avg. reaction rate in link
    FRICTION,      // link friction factor
  
    POWER,         // pump power output
    TIME,          // simulation time
    VOLUME,        // tank volume
    CLOCKTIME,     // simulation time of day
    FILLTIME,      // time to fill a tank
    DRAINTIME,     // time to drain a tank
    MAXVAR         // total number of variable fields
}

enum SectionType {
    _TITLE, _JUNCTIONS, _RESERVOIRS, _TANKS, _PIPES, _PUMPS,
    _VALVES, _CONTROLS, _RULES, _DEMANDS, _SOURCES, _EMITTERS,
    _PATTERNS, _CURVES, _QUALITY, _STATUS, _ROUGHNESS, _ENERGY,
    _REACTIONS, _MIXING, _REPORT, _TIMES, _OPTIONS,
    _COORDS, _VERTICES, _LABELS, _BACKDROP, _TAGS, _END
}

enum HdrType {
    STATHDR,       // hydraulic status header
    ENERHDR,       // energy usage header
    NODEHDR,       // node results header
    LINKHDR        // link results header
}

pub enum FlowDirection {
    NEGATIVE  = -1,  // flow in reverse of pre-assigned direction
    ZERO_FLOW = 0,   // zero flow
    POSITIVE  = 1    // flow in pre-assigned direction
}

enum DemandModelType {
    DDA,           // demand driven analysis
    PDA            // pressure driven analysis
}


/*
------------------------------------------------------
   Fundamental Data Structures
------------------------------------------------------
*/

pub struct IDString {
    ID: String
}

pub struct Spattern             // Time Pattern Object
{
  ID: String,             // pattern ID
  Comment: String,         // pattern comment
  Length: i32,          // pattern length
  F: Option<Vec<f64>>           // pattern factors
}

pub struct Scurve {
    ID: String,   // curve ID
    Comment: String,      // curve comment
    Type: CurveType,          // curve type
    Npts: i32,          // number of points
    Capacity: i32,      // size of X & Y arrays
    X: Option<Vec<f64>>,            // x-values
    Y: Option<Vec<f64>>            // y-values
}

pub struct Sdemand {
    Base: f64,             // baseline demand
    Pat: i32,              // pattern index
    Name: String,           // demand category name
    next: Option<Box<Sdemand>>    // next demand list item
}

pub struct Senergy {
    TimeOnLine: f64,       // hours pump is online
    Efficiency: f64,       // total time wtd. efficiency
    KwHrsPerFlow: f64,     // total kw-hrs per unit of flow
    KwHrs: f64,            // total kw-hrs consumed
    MaxKwatts: f64,        // max. kw consumed
    TotalCost: f64,        // total pumping cost
    CurrentPower: f64,     // current pump power (kw)
    CurrentEffic: f64     // current pump efficiency
}

pub struct Ssource {
    C0: f64,         // base concentration/mass
    Pat: i32,        // pattern index
    Smass: f64,      // actual mass flow rate
    Type: SourceType       // type of source
}

pub struct Svertices {
    X: Box<f64>,               // array of x-coordinates
    Y: Box<f64>,               // array of y-coordinates
    Npts: i32,             // number of vertex points
    Capacity: i32         // capacity of coordinate arrays
}

pub struct Snode {
    ID: String,    // node ID
    X: f64,              // x-coordinate
    Y: f64,              // y-coordinate
    El: f64,             // elevation
    D: Box<Sdemand>,              // demand pointer
    S: Box<Ssource>,              // source pointer
    C0: f64,             // initial quality
    Ke: f64,             // emitter coeff.
    Rpt: i32,            // reporting flag
    ResultIndex: i32,    // saved result index
    Type: NodeType,           // node type
    Comment: String       // node comment
}

pub struct Slink {
    ID: String,    // link ID
    N1: i32,             // start node index
    N2: i32,             // end node index
    Diam: f64,           // diameter
    Len: f64,            // length
    Kc: f64,             // roughness
    Km: f64,             // minor loss coeff.
    Kb: f64,             // bulk react. coeff.
    Kw: f64,             // wall react. coef.
    R: f64,              // flow resistance
    Rc: f64,             // reaction coeff.
    Type: LinkType,           // link type
    Status: StatusType,       // initial status
    VerticesL: Box<Svertices>,     // internal vertex coordinates
    Rpt: i32,            // reporting flag
    ResultIndex: i32,    // saved result index
    Comment: String       // link comment
}

pub struct Stank {
    Node: i32,            // node index of tank
    A: f64,               // tank area
    Hmin: f64,            // minimum water elev
    Hmax: f64,            // maximum water elev
    H0: f64,              // initial water elev
    Vmin: f64,            // minimum volume
    Vmax: f64,            // maximum volume
    V0: f64,              // initial volume
    Kb: f64,              // bulk reaction coeff.
    V: f64,               // tank volume
    C: f64,               // concentration
    Pat: i32,             // fixed grade time pattern
    Vcurve: i32,          // volume v. elev. curve index
    MixModel: MixType,        // type of mixing model
    V1frac: f64,          // mixing compartment fraction
    CanOverflow: i32     // tank can overflow or not
}

pub struct Spump {
    Link: i32,            // link index of pump
    Ptype: i32,           // pump curve type
    Q0: f64,              // initial flow
    Qmax: f64,            // maximum flow
    Hmax: f64,            // maximum head
    H0: f64,              // shutoff head
    R: f64,               // flow coeffic.
    N: f64,               // flow exponent
    Hcurve: i32,          // head v. flow curve index
    Ecurve: i32,          // effic. v. flow curve index
    Upat: i32,            // utilization pattern index
    Epat: i32,            // energy cost pattern index
    Ecost: f64,           // unit energy cost
    Energy: Senergy,          // energy usage statistics
}

pub struct Svalve {
    Link: i32 //  link index of valve
}

pub struct Scontrol {
    Link: i32,      // link index
    Node: i32,      // control node index
    Time: u32,      // control time
    Grade: f64,     // control grade
    Setting: f64,   // new link setting
    Status: StatusType,    // new link status
    Type: ControlType      // control type
}

pub struct SField {
    Name: String,  // name of reported variable
    Units: String, // units of reported variable
    Enabled: i32,        // enabled if in table
    Precision: i32,      // number of decimal places
    RptLim: [f64;2]      // lower/upper report limits
}

#[derive(Debug, PartialEq, Eq, Hash,)]
pub struct Sadjlist {
    pub node: i32,           // index of connecting node
    pub link: i32,           // index of connecting link
}

pub type Padjlist = Box<List<Sadjlist>>;

#[derive(Debug, PartialEq, Eq, Hash)]
pub struct Sseg {
    v: f64,             // segment volume
    c: f64,             // segment water quality
    // @FIXME: need prev instead of next, so 
    // List should be tweeked to provide the same
    // Capability
    // prev: Option<Box<Sseg>>    // previous segment in list
}

// @FIXME: Notice the Sseg for prev pointer instead of next
pub type Pseg = List<Sadjlist>;

#[derive(Debug, PartialEq, Eq, Hash)]
pub struct Spremise {
    logop: i32,            // logical operator (IF, AND, OR)
    object: i32,           // NODE or LINK
    index: i32,            // object's index
    variable: i32,         // pressure, flow, etc.
    relop: i32,            // relational operator (=, >, <, etc.)
    status: i32,           // variable's status (OPEN, CLOSED)
    value: f64,            // variable's value
}

#[derive(Debug, PartialEq, Eq, Hash)]
pub struct Saction {
    link: i32,              // link index
    status: i32,            // link's status
    setting: f64,           // link's setting
}

pub struct Srule {
    label: String,   // rule label
    priority: f64,         // priority level
    Premises: List<Spremise>,        // list of premises
    ThenActions: List<Saction>,     // list of THEN actions
    ElseActions: List<Saction>     // list of ELSE actions
}

#[derive(Debug, PartialEq, Eq, Hash,)]
pub struct SactionItem {
    ruleIndex: i32,           // index of rule action belongs to
    action: Option<Box<Saction>>,            // an action clause
}

pub type SactionList = List<SactionItem>

pub struct SmassBalance {
    initial: f64,         // initial mass in system
    inflow: f64,          // mass inflow to system
    outflow: f64,         // mass outflow from system
    reacted: f64,         // mass reacted in system
    r#final: f64,           // final mass in system
    ratio: f64,           // ratio of mass added to mass lost
}

/*
------------------------------------------------------
  Wrapper Data Structures
------------------------------------------------------
*/

// Input file parser wrapper
pub struct Parser {

    InFile: Box<File>,            // Input file handle

    DefPatID: String,     // Default demand pattern ID
    InpFname: String,  // Input file name
    Tok: Option<Vec<String>>,           // Array of token strings
    Comment: String,     // Comment text
    LineComment: String, // Full line comment

  
    MaxNodes: i32,              // Node count   from input file */
    MaxLinks: i32,              // Link count    "    "    "
    MaxJuncs: i32,              // Junction count "   "    "
    MaxPipes: i32,              // Pipe count    "    "    "
    MaxTanks: i32,              // Tank count    "    "    "
    MaxPumps: i32,              // Pump count    "    "    "
    MaxValves: i32,             // Valve count   "    "    "
    MaxControls: i32,           // Control count "   "     "
    MaxRules: i32,              // Rule count    "   "     "
    MaxPats: i32,               // Pattern count "   "     "
    MaxCurves: i32,             // Curve count   "   "     "
    Ntokens: i32,               // Number of tokens in line of input
    Ntitle: i32,                // Number of title lines
    ErrTok: i32,                // Index of error-producing token
    Unitsflag: i32,             // Unit system flag
    Flowflag: i32,              // Flow units flag
    Pressflag: i32,             // Pressure units flag
    DefPat: i32,                // Default demand pattern

    PrevPat: Option<Box<Spattern>>,       // Previous pattern processed
    PrevCurve: <Box<Scurve>>,     // Previous curve processed
    X: Option<Vec<f64>>,        // Temporary array for curve data
}

// Time step warpper
// @TODO: rust provide richer api for time management
// we might need to optimize this in future
pub struct Times {  
    Tstart: f32,                // Starting time of day
    Hstep: f32,                 // Nominal hyd. time step
    Pstep: f32,                 // Time pattern time step
    Pstart: f32,                // Starting pattern time
    Rstep: f32,                 // Reporting time step
    Rstart: f32,                // Time when reporting starts
    Rtime: f32,                 // Next reporting time
    Htime: f32,                 // Current hyd. time
    Hydstep: f32,               // Actual hydraulic time step
    Qstep: f32,                 // Quality time step
    Qtime: f32,                 // Current quality time
    Rulestep: f32,              // Rule evaluation time step
    Dur: f32                 // Duration of simulation
}

// Reporting warrper
pub struct Report {

    RptFile: Option<Box<File>>,           // Report file handle
    
    Nperiods:  i32,              // Number of reporting periods
    PageSize:  i32,              // Lines/page in output report/
    Rptflag:  i32,               // Report flag
    Tstatflag:  i32,             // Report time series statistic flag
    Summaryflag:  i32,           // Report summary flag
    Messageflag:  i32,           // Error/warning message flag
    Statflag:  i32,              // Status report flag
    Energyflag:  i32,            // Energy report flag
    Nodeflag:  i32,              // Node report flag
    Linkflag:  i32,              // Link report flag
    Fprinterr:  i32,             // File write error flag
    
    LineNum: f32,               // Current line number
    PageNum: f32,               // Current page number
    
    // @TODO: optimize this by using rust timestamp module
    Atime: [u8; 13],             // Clock time (hrs:min:sec)
    Rpt1Fname: String, // Primary report file name
    Rpt2Fname: String, // Secondary report file name
    DateStamp: [u8; 26],         // Current date & time
  
    Field: [SField; FieldType::MAXVAR as usize] // Output reporting fields
}

pub struct Outfile {
    
    HydFname: String,  // Hydraulics file name
    OutFname: String,  // Binary output file name

    Outflag: i32,               // Output file flag
    Hydflag: i32,               // Hydraulics flag
    SaveHflag: i32,             // Hydraulic results saved flag
    SaveQflag: i32,             // Quality results saved flag
    Saveflag: i32,              // General purpose save flag

  
    HydOffset: f32,             // Hydraulics file byte offset
    OutOffset1: f32,            // 1st output file byte offset
    OutOffset2: f32,            // 2nd output file byte offset

  
    OutFile: Option<Box<File>>,              // Output file handle
    HydFile: Option<Box<File>>,              // Hydraulics file handle
    TmpOutFile: Option<Box<File>>,         // Temporary file handle
}

// Rule-Based Controls Wrapper
pub struct Rules {
    
    ActionList: Option<Box<SactionList>>,     // Linked list of action items
    RuleState: i32,       // State of rule interpreter
    Errcode: i32,         // Rule parser error code
    Time1: f32,           // Start of rule evaluation time interval
    LastPremise: Option<Box<Spremise>>,    // Previous premise clause
    LastThenAction: Option<Box<Saction>>, // Previous THEN action
    LastElseAction: Option<Box<Saction>>, // Previous ELSE action
}

// Sparse matrix wrapper
pub struct Smatrix {
    
    pub Aii: Option<Vec<f32>>,        // Diagonal matrix coeffs.
    pub Aij: Option<Vec<f32>>,        // Non-zero, off-diagonal matrix coeffs.
    pub F: Option<Vec<f32>>,          // Right hand side vector
    pub temp: Option<Vec<f32>>,       // Array used by linear eqn. solver

    pub Ncoeffs: i32,    // Number of non-zero matrix coeffs
    pub Order: Option<Vec<i32>>,      // Node-to-row of re-ordered matrix
    pub Row: Option<Vec<i32>>,        // Row-to-node of re-ordered matrix
    pub Ndx: Option<Vec<i32>>,        // Index of link's coeff. in Aij
    pub XLNZ: Option<Vec<i32>>,       // Start position of each column in NZSUB
    pub NZSUB: Option<Vec<i32>>,      // Row index of each coeff. in each column
    pub LNZ: Option<Vec<i32>>,        // Position of each coeff. in Aij array
    pub link: Option<Vec<i32>>,       // Array used by linear eqn. solver
    pub first: Option<Vec<i32>>,      // Array used by linear eqn. solver
}

// Hydraulics solver wrapper
pub struct Hydraul {
    
    NodeHead: Option<Vec<f32>>,             // Node hydraulic heads
    NodeDemand: Option<Vec<f32>>,           // Node demand + emitter flows
    DemandFlow: Option<Vec<f32>>,           // Work array of demand flows
    EmitterFlow: Option<Vec<f32>>,          // Emitter outflows
    LinkFlow: Option<Vec<f32>>,             // Link flows
    LinkSetting: Option<Vec<f32>>,          // Link settings
    Htol: f32,                  // Hydraulic head tolerance
    Qtol: f32,                  // Flow rate tolerance
    RQtol: f32,                 // Flow resistance tolerance
    Hexp: f32,                  // Exponent in headloss formula
    Qexp: f32,                  // Exponent in emitter formula
    Pexp: f32,                  // Exponent in demand formula
    Pmin: f32,                  // Pressure needed for any demand
    Preq: f32,                  // Pressure needed for full demand
    Dmult: f32,                 // Demand multiplier
    Hacc: f32,                  // Relative flow change limit
    FlowChangeLimit: f32,       // Absolute flow change limit
    HeadErrorLimit: f32,        // Hydraulic head error limit
    DampLimit: f32,             // Solution damping threshold
    Viscos: f32,                // Kin. viscosity (sq ft/sec)
    SpGrav: f32,                // Specific gravity
    Epump: f32,                 // Global pump efficiency
    Dsystem: f32,               // Total system demand
    Ecost: f32,                 // Base energy cost per kwh
    Dcost: f32,                 // Energy demand charge/kw/day
    Emax: f32,                  // Peak energy usage
    RelativeError: f32,         // Total flow change / total flow
    MaxHeadError: f32,          // Max. error for link head loss
    MaxFlowChange: f32,         // Max. change in link flow
    DemandReduction: f32,       // % demand reduction at pressure deficient nodes
    RelaxFactor: f32,           // Relaxation factor for flow updating
    P: Box<f32>,                    // Inverse of head loss derivatives
    Y: Box<f32>,                    // Flow correction factors
    Xflow: Box<f32>,                // Inflow - outflow at each node

    Epat: i32,                  // Energy cost time pattern
    DemandModel: i32,           // Fixed or pressure dependent
    Formflag: i32,              // Head loss formula flag
    Iterations: i32,            // Number of hydraulic trials taken
    MaxIter: i32,               // Max. hydraulic trials allowed
    ExtraIter: i32,             // Extra hydraulic trials
    CheckFreq: i32,             // Hydraulic trials between status checks
    MaxCheck: i32,              // Hydraulic trials limit on status checks
    OpenHflag: i32,             // Hydraulic system opened flag
    Haltflag: i32,              // Flag to halt simulation
    DeficientNodes: i32,        // Number of pressure deficient nodes

  
    LinkStatus: Option<Vec<StatusType>>,           // Link status
    OldStatus: Option<Vec<StatusType>>,            // Previous link/tank status
    
    pub smatrix: Smatrix         // Sparse matrix storage

}

// Water quality solver wrapper
pub struct Quality {

  Qualflag: i32,              // Water quality analysis flag
  OpenQflag: i32,             // Quality system opened flag
  Reactflag: i32,             // Reaction indicator
  OutOfMemory: i32,           // Out of memory indicator
  TraceNode: i32,             // Source node for flow tracing
  SortedNodes: Box<i32>,         // Topologically sorted node indexes

  pub ChemName : String,   // Name of chemical
  pub ChemUnits : String,  // Units of chemical

  Ctol: f64,                  // Water quality tolerance
  Diffus: f64,                // Diffusivity (sq ft/sec)
  Wbulk: f64,                 // Avg. bulk reaction rate
  Wwall: f64,                 // Avg. wall reaction rate
  Wtank: f64,                 // Avg. tank reaction rate
  Wsource: f64,               // Avg. mass inflow
  Rfactor: f64,               // Roughness-reaction factor
  Sc: f64,                    // Schmidt Number
  Bucf: f64,                  // Bulk reaction units conversion factor
  Tucf: f64,                  // Tank reaction units conversion factor
  BulkOrder: f64,             // Bulk flow reaction order
  WallOrder: f64,             // Pipe wall reaction order
  TankOrder: f64,             // Tank reaction order
  Kbulk: f64,                 // Global bulk reaction coeff.
  Kwall: f64,                 // Global wall reaction coeff.
  Climit: f64,                // Limiting potential quality
  SourceQual: f64,            // External source quality
  NodeQual: Option<Vec<f64>>,             // Reported node quality state
  PipeRateCoeff: Option<Vec<f64>>,        // Pipe reaction rate coeffs.

  FirstSeg: Option<Box<Box<Sseg>>>,             // First (downstream) segment in each pipe
  LastSeg: Option<Box<Sseg>>,              // Last (upstream) segment in each pipe
  FlowDir: Option<Box<FlowDirection>>,              // Flow direction for each pipe

  MassBalance: SmassBalance         // Mass balance components
}

// Pipe entwork wrapper
pub struct Network {
    
    pub Nnodes: i32,                // Number of network nodes
    Ntanks: i32,                // Number of tanks
    Njuncs: i32,                // Number of junction nodes
    pub Nlinks: i32,                // Number of network links
    Npipes: i32,                // Number of pipes
    Npumps: i32,                // Number of pumps
    Nvalves: i32,               // Number of valves
    Ncontrols: i32,             // Number of simple controls
    Nrules: i32,                // Number of control rules
    Npats: i32,                 // Number of time patterns
    Ncurves: i32,               // Number of data curves

    Node: Option<Vec<Snode>>,          // Node array
    Link: Option<Vec<Slink>>,          // Link array
    Tank: Option<Vec<Stank>>,          // Tank array
    Pump: Option<Vec<Spump>>,          // Pump array
    Valve: Option<Vec<Svalve>>,         // Valve array
    Pattern: Option<Vec<Spattern>>,       // Time pattern array
    Curve: Option<Vec<Scurve>>,         // Data curve array
    Control: Option<Vec<Scontrol>>,       // Simple controls array
    Rule: Option<Vec<Srule>>,          // Rule-based controls array
    
    NodesHashMapTable: HashMap<String, String>,        // Hash Map for Node ID names
    LinkHashTable: HashMap<String, String>,        // Hash Map for Link ID names
    // @FIXME: should make sure that this
    pub Adjlist: Option<Box<Padjlist>>       // Node adjacency lists
}


pub struct Project {
    pub network: Network,            // Pipe network wrapper
    pub parser: Parser ,             // Input file parser wrapper
    times: Times  ,              // Time step wrapper
    report: Report ,             // Reporting wrapper
    outfile: Outfile,            // Output file wrapper
    rules: Rules  ,              // Rule-based controls wrapper
    pub hydraul: Hydraul,            // Hydraulics solver wrapper
    quality: Quality,            // Water quality solver wrapper

    Ucf: [f64; FieldType::MAXVAR as usize],            // Unit conversion factors

    Openflag: i32,                    // Project open flag
    Warnflag: i32,                    // Warning flag

    Msg: String,             // General-purpose string: errors, messages
    Title: String,           // Project title
    MapFname: String,        // Map file name
    TmpHydFname: String,     // Temporary hydraulics file name
    TmpOutFname: String,     // Temporary output file name
    TmpStatFname: String     // Temporary statistic file name
}
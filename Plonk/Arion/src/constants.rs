use crate::{WIDTH, ROUNDS, MAX_WIDTH, MAX_ROUNDS};
use dusk_bls12_381::BlsScalar;

pub const MATRIX: [[BlsScalar; WIDTH]; WIDTH] = {
    let mut matrix = [[BlsScalar::zero(); WIDTH]; WIDTH];

    let one = BlsScalar::one();

    let mut i = 1;
    while i < WIDTH {
        matrix[0][i] = matrix[0][i  - 1].add(&one);
        i += 1;
    }

    i = 1;
    let mut j = 1;
    while i < WIDTH {
        matrix[i][0] = matrix[i - 1][WIDTH - 1];
        while j < WIDTH {
            matrix[i][j] = matrix[i - 1][j - 1];
            j += 1;
        }
        j = 1;
        i += 1;
    }

    matrix
};

pub const CONSTANTS_AFF: [[BlsScalar; WIDTH]; ROUNDS] = {
    let const_array: [[[u64; 4]; MAX_WIDTH]; MAX_ROUNDS] = [
    [[970454385332725214, 4735158914852834846, 16702983313962767143, 7375736764391876140], [4230679363595362084, 14748170609639716385, 6726276718525653011, 209585127932669456], [6985368126114302085, 15042707609042090796, 6574303427188805453, 4260491623787833808], [4211534376445758016, 4194194800130890524, 6189979315052828417, 6333043189070472884], [7713061483966079236, 2531813565031281972, 15141872849730284883, 2617614298691088635], [9384435332689772208, 6380797830192786222, 14046751767082687730, 6726703997800057086], [1992080556398414636, 4864203202662323702, 12166009503139942213, 4341074411686198705], [9183113592048845618, 18038533873937178850, 16188992193114910863, 5450448736873179721]],
    [[1266687978422234986, 69881920587063368, 16171726449901498918, 5950717589371861074], [6930917184652543673, 10514231817671922317, 4035945929007950627, 6501921180983397386], [7732435664247685777, 1843882157446206982, 5208046253713768039, 3231264108358203476], [2468010120310035552, 752360851117569350, 16539927831095994414, 4965934128765020503], [13585508233168112285, 14268525048955505113, 3893527268571935914, 6864401681748406195], [7757576072757344606, 5938280521025358971, 13307080065139554079, 7861671533100997204], [16975249928919073494, 14924052493920374330, 8467595251003959199, 4341880216694556651], [1632378518933054389, 5376691641771404223, 6426806902617696306, 3820596342789430525]],
    [[3554599183847677264, 6791278616554961889, 63130878131432096, 4185140968160079833], [8144674920183683498, 8216212977195016227, 9686337916691652703, 6076711101860273494], [12922228869988892970, 15727682283054545143, 2485208998993290360, 6509977704032023471], [7539704408897155719, 18446043554162454165, 4910234100827969386, 7967271706650199014], [2941992049881786276, 3273283349677813825, 9270159677016703884, 5881989563811884977], [15887736589393534182, 8055049948876732053, 966079056558760390, 6559689159628549478], [5353974172002425974, 2257867930741806795, 17066911527045229001, 5886847415621009780], [2689410499930379792, 13110773802046837508, 5774031095598821402, 1662368076039383236]],
    [[3571850794035492638, 10011780587604914872, 2851532754956463661, 7549769893150402939], [11021022793020711086, 3190937485316022483, 16963355507106829895, 4540483237117649868], [7950125929806482344, 4629806249716173727, 17762446330104727293, 4599927447489006632], [1973144313235787763, 4063288474533805473, 17925690664622908455, 1719901659978698241], [14191362668716585383, 5514923056418079799, 17106225004840810921, 1099833732327433935], [17249587290320992629, 15842981470461115274, 6532834680074548848, 1787665309043674696], [14183316565875525366, 5704389414999837909, 4997229107221395814, 8059036405762643152], [14342287689315652981, 13433382801171108475, 3238436758015374989, 6216378968260281689]],
    [[9013631226858709030, 14343543889092537838, 17937032940367356980, 5584102768109416368], [9097166776537815687, 11222723496119898897, 1023023202605269891, 7766214631594472005], [9073964238643267673, 12299973005925769039, 4683437631921202259, 6944843882550015661], [15162245261304966082, 7473032242873636412, 8110499706316289816, 818088745158106631], [6384609283665166189, 5213498141675544488, 16564560986343493811, 7130572796089109734], [2614091791817208963, 1713030340148174225, 11222282542926130617, 1252455153930172006], [3927563877318024503, 12092270890523652051, 7276647802898046790, 7025383794357204764], [4894603339441868829, 16094733858056055320, 4423825234437775665, 3126972067811597095]],
    [[10147312101186601355, 16142173913991205571, 9954551831634786764, 499905868224188646], [12820247075599576363, 3149924745696910946, 4963691744052408888, 5868682240901226497], [15827773287667575720, 14213358503609674465, 10560250937334462856, 4386318857047935211], [9397446180221503578, 6618786245379209115, 3084246448925900426, 2822857807936729361], [429493542051960951, 6682543568550011335, 14477988022801782791, 6477944054337417900], [13530363729867616788, 3022258258948725978, 17074046852931477145, 6762314476481292508], [15289080667886528133, 13249111072803844268, 13125954514015266183, 5252959941099154400], [7426130465162481050, 909871896371706739, 26466638563414140, 7269823781459537969]],
    ];

    let mut affine_constants = [[BlsScalar::zero(); WIDTH]; ROUNDS];

    let mut i = 0;
    let mut j = 0;
    while i < ROUNDS {
        while j < WIDTH {
            affine_constants[i][j] = BlsScalar::from_raw(const_array[i][j]);
            j += 1;
        }
        j = 0;
        i += 1;
    }

    affine_constants
};

pub const CONSTANTS_G: [[[BlsScalar; 2]; WIDTH - 1]; ROUNDS] = {
    let const_array: [[[[u64; 4]; 2]; MAX_WIDTH - 1]; MAX_ROUNDS] = [
        [[[8118021754834449418, 18137426088772698758, 11203759450720152331, 3900225123309009714], [14639014844555717986, 10036122947945021134, 16847282596167317246, 7187814803256400694]], [[14525900076367906240, 285701828013040820, 8209039863985472825, 6002201468642717514], [1976020842633079501, 10845738633802611111, 13765667994085148964, 3261502654422324242]], [[7758101740920843174, 18176659492512981487, 8932647722209367966, 6536481862122956546], [18144900633143714544, 9006025712098685260, 17305636935215169546, 4021291266829413815]], [[13430340930673190471, 6340972178702447167, 18395728821334867436, 4507679155598365721], [13274273673115050468, 14355740628689897277, 16358078499105545713, 4206170965530390973]], [[12969467437587247133, 14185495908945610644, 3356607114650797197, 6562861671434439], [5288709648859117752, 7862515103577432033, 17470788809537922822, 5858991881344592831]], [[15909451001460723575, 15899331260159970773, 12497740786161253656, 857273899484426815], [1726323722231622433, 14164303893094005697, 2767468258502354084, 3639240972191378544]], [[10685093130424192313, 4650403758870880299, 15817251630106561038, 6819670475876146587], [15482540654177638685, 8173073208326405382, 8432768059150967505, 2329304385750892884]]],
        [[[9545677499604737483, 9462316435496250230, 6219837032005400383, 1669594911932636079], [14667649992759821223, 15250828530975944016, 12461680458907490591, 2600250255282402069]], [[13479747965968902564, 11521793992377043895, 16908475898601076588, 5922391494409083744], [12830384399261354384, 7405181226035967048, 7340916367803008605, 5212575767673113425]], [[5937793381302090817, 7415845228460016657, 17925875622089189925, 6165785024554652084], [12506554955952508363, 9615210230693738864, 4754073028702950625, 4451246544233010362]], [[15486209265845486165, 7654925983749355124, 13473162313556685240, 1997024781483169030], [17581342995316276378, 1113040793097221465, 3686486519207174303, 7415607786627847148]], [[16965840722640656121, 14525749725805698566, 544683941451307980, 4425297321159290825], [8291601896729146924, 4779904461668358139, 16247936753260975579, 3039784398180502571]], [[12047240076415696150, 13759167846029586356, 6399470930210199921, 1714357474642878917], [15108671714317370091, 10194000484015189560, 11177671332330699425, 3737361221871765393]], [[3060813626616936872, 7959071218933046848, 7662273436732224696, 5309218351774726885], [2849143678355831605, 3574989035684888073, 1175389452297139059, 8108499115226511546]]],
        [[[208280022475034388, 2679639529983808676, 12013685850765403359, 3129571291788010104], [15033347828443971160, 16172786514232009685, 5296479324246923749, 7916974900529105519]], [[6403808021927291531, 1882243198164585761, 10747951411459338398, 5286309658413882805], [2589338544038484845, 892230472111265022, 9176940993904317953, 7900007709313640766]], [[5825135915400858962, 16456310909329014853, 3688649171574841328, 6937065908095652933], [7921322509731799822, 5496356022891403552, 428584605920154781, 5132614768929143708]], [[15067970489178591057, 600730711982177270, 1590253018314966921, 968851811485403290], [16547425915152315840, 10532392093557729319, 3404982871926017629, 5740444085366615395]], [[17917962971050611055, 4895448124054798493, 5753427114127888753, 743770653188465394], [6943716825838383689, 14171567719930357365, 10250355854456927139, 6363648223488998798]], [[10459870108950237175, 1693082505502488367, 14945752005013406473, 1366982426230626313], [6532539179302865979, 12425412972709486351, 13197308402304315423, 3799031331348979717]], [[1380350282873968945, 5395281819162148170, 6552630896469362116, 2358421546768082294], [14989364060735952308, 891855254144811967, 16194507144776794135, 4264912691036323915]]],
        [[[6353085445362285923, 7658188572571943832, 6574250057995045111, 7768576349785385851], [7648284924978271806, 2877333109811476230, 6498093629955401021, 7019283868739783883]], [[17056219359150660356, 5196123004898719298, 12375888835881030296, 3612529847840289789], [1076413377290184717, 18360863644199700368, 6964199429652363574, 4181798701167475269]], [[14314825351960927816, 9919977589572721358, 9252906103615682771, 6251822302271510673], [4250568383850889233, 4926391120193127295, 16540102346049935944, 7948619022851418508]], [[13302572416082948905, 4909947022186490706, 5895484041106108370, 6692524739130276136], [7009425706915053159, 17988128018894260530, 11052494136951769444, 6858228907836942337]], [[16030470096622617676, 12568212927151885586, 6186168307824137253, 7111668131997318839], [6390845260264110765, 4126870813794157858, 3059978959119024563, 4576934722927331636]], [[16239501893001237546, 7059152550827108346, 3685693859866159845, 3478486378672058961], [10882933867352097995, 9357076616662967518, 14026968747970949809, 4939422518531225136]], [[2718922157413486096, 9110711474542097594, 7637475741053790091, 1957593766721549182], [5443882364869207071, 11670732637609796819, 1997461725720262185, 7984598109576538120]]],
        [[[12406281688707078852, 9395863153120753303, 3299817659868299112, 3307149134112359590], [7828540874485720701, 13486290082861638100, 11281604080227608769, 496391192300390777]], [[7862386140070849403, 15289371184814986741, 14049346157782498755, 1993375538878569720], [10764251912436354296, 15272485557701176198, 2551461113230616347, 4514633682531610571]], [[12804934626368819242, 7739524292090359383, 3468334448200733679, 6322067415189582302], [2229713910601436735, 6345558670892006099, 17205148658989621795, 3884981506912411636]], [[12167776622922715529, 14734532337811862294, 6330850869321644243, 2491935068974958207], [11170561174960155587, 8117455006381980714, 6808988857912942901, 1354717038749667982]], [[3525528504342493763, 8928014825702114870, 119202813826337560, 7348215879624833139], [13948477937432111501, 82835732227396954, 12908772836769617377, 3618434982189587858]], [[12452690745000238041, 4014262368108768616, 13857388782199059524, 7354634678231609479], [11604181008145272736, 13278800694907059260, 945587509935036877, 7132176450243810895]], [[5275625478783838709, 7529966767911732600, 6173886294747637985, 112760325689742107], [9844793801479437508, 17152235997821162554, 994191905673958105, 5517555803449611114]]],
        [[[12190053303382890390, 16700011988895785313, 15144368542008125936, 238249394300578986], [8979763215614725737, 15585637802562290449, 11707763093106247433, 1453886827838719509]], [[15958231742389272731, 3320088711834700778, 14508091930884470503, 5996888392890636072], [17558997193776029659, 8517004962181611079, 18323119281018825345, 4072424636002956544]], [[8224770923119548965, 14987223647211525204, 9512075764922500944, 5379655282320561415], [2121680431317067555, 5227943416316865675, 15086880451186976961, 647184733967484958]], [[9171366037063595778, 3798752676398428079, 13275376650134637232, 8306416127283865314], [17267067399339640886, 15458458084302052279, 12014148314634911025, 5259416227506602583]], [[2319639134903849842, 6530716520768153299, 10496585537729916184, 6993254964715602817], [14865093722343191197, 17233573779832860554, 6918475126755900122, 208461796375197346]], [[14862115151621272843, 9626950628945589338, 15666928646406271507, 6973139649509326923], [4079848560759088650, 345977733172885035, 17120375523020352721, 4972476976776665811]], [[2757805764682775139, 3122048073449868516, 3752364280538773500, 5729254583848923047], [7738468962757172082, 9631403622297223643, 17966057605146409821, 8325318453324702420]]],
        ];

    let mut constants_g = [[[BlsScalar::zero(); 2]; WIDTH - 1]; ROUNDS];

    let mut i = 0;
    let mut j = 0;
    let mut tmp = [BlsScalar::zero(); 2];
    while i < ROUNDS {
        while j < WIDTH - 1{
            tmp[0] = BlsScalar::from_raw(const_array[i][j][0]);
            tmp[1] = BlsScalar::from_raw(const_array[i][j][1]);
            constants_g[i][j] = tmp;
            j += 1;
        }
        j = 0;
        i += 1;
    }

    constants_g
};

pub const CONSTANTS_H: [[BlsScalar; WIDTH - 1]; ROUNDS] = {
    let const_array: [[[u64; 4]; MAX_WIDTH - 1]; MAX_ROUNDS] = [
        [[10921243481068140754, 2518525656384150226, 6113210645428543405, 3965711976898692104], [3485595414554252575, 409582036943307088, 13828756501544561109, 1677462619007352731], [3334782092255076761, 5570395682066137536, 10608845724528439175, 7579215251462115483], [9694456976992470658, 15804342176037893448, 9169995828282038490, 3088629552174663926], [4736994860782243944, 7685556290627959999, 11715635745556264818, 5943033524509355117], [10215347721694543935, 15993796356778898738, 10137559064734388862, 227314686331948528], [13460852012916779609, 4358098711408330618, 8978707826599623711, 7213924270875377657]],
        [[4926870207870077544, 5277045161118171188, 10224442989354951543, 7444723940524320645], [7471186748802007366, 7260661708697062126, 7283265288265511831, 1869305939512461199], [3095833256795908432, 16687801437939668899, 7117743567829038038, 3849376791502321111], [7960370455870505448, 16454464364719604501, 4344917559538738003, 8090736631896107454], [13146371958286040330, 14077105886872986004, 15014754281816200845, 8026026921204648236], [10302803153185923433, 3785189656327363139, 4070816978900273521, 6377282511292992837], [4568675908884195690, 5719935760946173319, 7853822667326590627, 1046559994850697583]],
        [[9319126307559398719, 7996037948054294134, 16392394823352622401, 8188308676097467852], [8039360918536705859, 8938748092297590930, 10099583915951756141, 952702014117514948], [10285667626563052611, 12135177443099297491, 17994610563755969311, 340137640551325376], [12510383504830329153, 8962682242325219484, 18068237144866482128, 6009124256220310615], [4979681699100286651, 13057368319490099009, 9302418961475703129, 2305390495217965735], [9880244046737069480, 15005119545105405351, 10135635245795893073, 2343835903853041702], [4720459097565055526, 41566771599115942, 5431896585204844030, 900462956678453416]],
        [[8180551441939956250, 1372082306606461051, 10939553004169839157, 4174004388515504730], [10781897078796356039, 11161264076082777382, 5358217306697285381, 1279740062567306272], [12880572968795548036, 778044184158216145, 7847837927596014841, 5093713696678773510], [10313630131122508629, 11595433834679230427, 17906853853013777625, 4485645081817223593], [13106500706271127444, 12887684643604891104, 8468942193355381622, 5996179271435188893], [11377160197169944980, 14200869053881868093, 6982204977217433805, 1071277133691504131], [7604845471956914678, 3109004676757811297, 7225943085490101937, 7919172033569140650]],
        [[7101833186995242907, 10107917178379272097, 17732701405176786982, 4154920907981975485], [1219830418552994049, 1369037724583717308, 2167715856742450452, 4631957203260248689], [10259558657442197957, 2192041270483450352, 14746560274696040402, 5119600357043938969], [10476521196996719913, 7784740041053923661, 6154106695954335971, 6413629037280236525], [17580664618409032270, 13690165057445147271, 5789944169470755250, 3769994094034477273], [12011514537309988796, 13631591614792306181, 15022188569403805110, 7451264850468272863], [5444580192748715976, 16484411435010131438, 658731927658763468, 1368787473530077100]],
        [[44232278147836182, 3395565509036815150, 10349233983607859680, 6182181519168384223], [17686354943380459589, 1869381740183272824, 11734889076207267566, 4351731507497101934], [7827189883797829130, 521796066750023363, 8782522078991107946, 1790836173019697558], [2422273977950441257, 7068580442719219078, 14611251581445613343, 3559553411702362473], [5260767482580517521, 8841526923719449880, 6329734848449174868, 6790249465063915528], [18409164497553330764, 17230383061068002211, 1748711811956533389, 5164881649052590803], [5674536440827512799, 4430291039549499721, 13334935409224003304, 7834257547345749173]],
    ];

    let mut constants_h = [[BlsScalar::zero(); WIDTH - 1]; ROUNDS];
    let mut i = 0;
    let mut j = 0;
    while i < ROUNDS {
        while j < WIDTH - 1{
            constants_h[i][j] = BlsScalar::from_raw(const_array[i][j]);
            j += 1;
        }
        j = 0;
        i += 1;
    }

    constants_h
};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{WIDTH, ROUNDS};

    fn legendre_symbol(val: BlsScalar) -> i8 {
        // (p - 1) / 2
        let exponent: [u64; 4] = [
            9223372034707292160,
            12240451741123816959,
            1845609449319885826,
            4176758429732224676,
        ];
        let one = BlsScalar::one();
        let zero = BlsScalar::zero();
        let symb = val.pow_vartime(&exponent);

        if symb == zero {
            0
        }
        else if symb == one {
            1
        }
        else {
            -1
        }
    }

    #[test]
    fn check_for_irreducibility_of_polynomials() {
        for i in 0..ROUNDS {
            for j in 0..(WIDTH - 1) {
                let val = CONSTANTS_G[i][j][0].square() - CONSTANTS_G[i][j][1].double().double();
                assert_eq!(legendre_symbol(val), -1 as i8);
            }
        }

    }
}
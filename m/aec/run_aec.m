close all;
clear all;
filter_coeff=[-0.000030328878, -0.000030577072, -0.000031071936, -0.000031810440, ...
        -0.000032788081, -0.000033998941, -0.000035435750, -0.000037089977, ...
        -0.000038951924, -0.000041010841, -0.000043255049, -0.000045672075, ...
        -0.000048248791, -0.000050971567, -0.000053826416, -0.000056799157, ...
        -0.000059875563, -0.000063041519, -0.000066283168, -0.000069587057, ...
        -0.000072940270, -0.000076330559, -0.000079746457, -0.000083177384, ...
        -0.000086613733, -0.000090046951, -0.000093469597, -0.000096875386, ...
        -0.000100259219, -0.000103617193, -0.000106946596, -0.000110245887, ...
        -0.000113514662, -0.000116753598, -0.000119964392, -0.000123149681, ...
        -0.000126312954, -0.000129458453, -0.000132591065, -0.000135716209, ...
        -0.000138839714, -0.000141967702, -0.000145106460, -0.000148262318, ...
        -0.000151441530, -0.000154650156, -0.000157893954, -0.000161178271, ...
        -0.000164507953, -0.000167887260, -0.000171319789, -0.000174808415, ...
        -0.000178355242, -0.000181961566, -0.000185627852, -0.000189353730, ...
        -0.000193137994, -0.000196978626, -0.000200872823, -0.000204817042, ...
        -0.000208807055, -0.000212838013, -0.000216904517, -0.000221000704, ...
        -0.000225120328, -0.000229256857, -0.000233403563, -0.000237553622, ...
        -0.000241700207, -0.000245836584, -0.000249956200, -0.000254052771, ...
        -0.000258120361, -0.000262153450, -0.000266147004, -0.000270096521, ...
        -0.000273998076, -0.000277848354, -0.000281644665, -0.000285384955, ...
        -0.000289067796, -0.000292692375, -0.000296258456, -0.000299766348, ...
        -0.000303216849, -0.000306611182, -0.000309950933, -0.000313237965, ...
        -0.000316474338, -0.000319662220, -0.000322803792, -0.000325901155, ...
        -0.000328956235, -0.000331970689, -0.000334945811, -0.000337882450, ...
        -0.000340780925, -0.000343640950, -0.000346461575, -0.000349241120, ...
        -0.000351977139, -0.000354666377, -0.000357304758, -0.000359887364, ...
        -0.000362408450, -0.000364861452, -0.000367239018, -0.000369533056, ...
        -0.000371734777, -0.000373834773, -0.000375823080, -0.000377689276, ...
        -0.000379422563, -0.000381011877, -0.000382445990, -0.000383713619, ...
        -0.000384803542, -0.000385704709, -0.000386406356, -0.000386898114, ...
        -0.000387170118, -0.000387213101, -0.000387018495, -0.000386578507, ...
        -0.000385886201, -0.000384935557, -0.000383721524, -0.000382240056, ...
        -0.000380488141, -0.000378463812, -0.000376166139, -0.000373595219, ...
        -0.000370752142, -0.000367638943, -0.000364258554, -0.000360614724, ...
        -0.000356711946, -0.000352555363, -0.000348150670, -0.000343504007, ...
        -0.000338621843, -0.000333510867, -0.000328177859, -0.000322629576, ...
        -0.000316872631, -0.000310913376, -0.000304757787, -0.000298411360, ...
        -0.000291879010, -0.000285164980, -0.000278272762, -0.000271205024, ...
        -0.000263963554, -0.000256549216, -0.000248961920, -0.000241200601, ...
        -0.000233263221, -0.000225146777, -0.000216847326, -0.000208360028, ...
        -0.000199679191, -0.000190798343, -0.000181710300, -0.000172407255, ...
        -0.000162880871, -0.000153122383, -0.000143122702, -0.000132872532, ...
        -0.000122362479, -0.000111583170, -0.000100525363, -0.000089180064, ...
        -0.000077538630, -0.000065592878, -0.000053335171, -0.000040758515, ...
        -0.000027856630, -0.000014624018, -0.000001056021, 0.000012851138, ...
        0.000027100327, 0.000041693485, 0.000056631617, 0.000071914809, ...
        0.000087542246, 0.000103512249, 0.000119822322, 0.000136469211, ...
        0.000153448970, 0.000170757041, 0.000188388336, 0.000206337332, ...
        0.000224598164, 0.000243164726, 0.000262030769, 0.000281190009, ...
        0.000300636218, 0.000320363327, 0.000340365516, 0.000360637298, ...
        0.000381173598, 0.000401969827, 0.000423021935, 0.000444326464, ...
        0.000465880589, 0.000487682136, 0.000509729602, 0.000532022152, ...
        0.000554559610, 0.000577342430, 0.000600371663, 0.000623648909, ...
        0.000647176257, 0.000670956213, 0.000694991628, 0.000719285612, ...
        0.000743841437, 0.000768662450, 0.000793751965, 0.000819113166, ...
        0.000844749009, 0.000870662116, 0.000896854682, 0.000923328383, ...
        0.000950084291, 0.000977122793, 0.001004443523, 0.001032045304, ...
        0.001059926097, 0.001088082964, 0.001116512045, 0.001145208547, ...
        0.001174166740, 0.001203379980, 0.001232840729, 0.001262540599, ...
        0.001292470407, 0.001322620233, 0.001352979500, 0.001383537057, ...
        0.001414281270, 0.001445200126, 0.001476281334, 0.001507512436, ...
        0.001538880923, 0.001570374340, 0.001601980408, 0.001633687126, ...
        0.001665482879, 0.001697356541, 0.001729297566, 0.001761296070, ...
        0.001793342912, 0.001825429752, 0.001857549108, 0.001889694393, ...
        0.001921859941, 0.001954041023, 0.001986233843, 0.002018435527, ...
        0.002050644089, 0.002082858396, 0.002115078108, 0.002147303617, ...
        0.002179535967, 0.002211776768, 0.002244028101, 0.002276292416, ...
        0.002308572425, 0.002340870987, 0.002373190996, 0.002405535266, ...
        0.002437906415, 0.002470306756, 0.002502738185, 0.002535202083, ...
        0.002567699223, 0.002600229679, 0.002632792754, 0.002665386912, ...
        0.002698009728, 0.002730657846, 0.002763326951, 0.002796011756, ...
        0.002828706000, 0.002861402467, 0.002894093007, 0.002926768579, ...
        0.002959419304, 0.002992034530, 0.003024602902, 0.003057112453, ...
        0.003089550693, 0.003121904707, 0.003154161264, 0.003186306922, ...
        0.003218328139, 0.003250211388, 0.003281943260, 0.003313510581, ...
        0.003344900504, 0.003376100618, 0.003407099027, 0.003437884440, ...
        0.003468446237, 0.003498774536, 0.003528860235, 0.003558695057, ...
        0.003588271569, 0.003617583196, 0.003646624221, 0.003675389765, ...
        0.003703875769, 0.003732078943, 0.003759996728, 0.003787627222, ...
        0.003814969120, 0.003842021626, 0.003868784372, 0.003895257323, ...
        0.003921440679, 0.003947334780, 0.003972940002, 0.003998256657, ...
        0.004023284898, 0.004048024621, 0.004072475378, 0.004096636296, ...
        0.004120505997, 0.004144082539, 0.004167363355, 0.004190345212, ...
        0.004213024178, 0.004235395601, 0.004257454099, 0.004279193573, ...
        0.004300607214, 0.004321687545, 0.004342426454, 0.004362815256, ...
        0.004382844754, 0.004402505313, 0.004421786948, 0.004440679412, ...
        0.004459172292, 0.004477255113, 0.004494917441, 0.004512148992, ...
        0.004528939733, 0.004545279991, 0.004561160551, 0.004576572754, ...
        0.004591508584, 0.004605960752, 0.004619922770, 0.004633389009, ...
        0.004646354755, 0.004658816245, 0.004670770698, 0.004682216325, ...
        0.004693152332, 0.004703578904, 0.004713497184, 0.004722909228, ...
        0.004731817958, 0.004740227096, 0.004748141095, 0.004755565048, ...
        0.004762504603, 0.004768965859, 0.004774955265, 0.004780479506, ...
        0.004785545399, 0.004790159771, 0.004794329357, 0.004798060681, ...
        0.004801359958, 0.004804232990, 0.004806685075, 0.004808720919, ...
        0.004810344567, 0.004811559332, 0.004812367747, 0.004812771523, ...
        0.004812771523, 0.004812367747, 0.004811559332, 0.004810344567, ...
        0.004808720919, 0.004806685075, 0.004804232990, 0.004801359958, ...
        0.004798060681, 0.004794329357, 0.004790159771, 0.004785545399, ...
        0.004780479506, 0.004774955265, 0.004768965859, 0.004762504603, ...
        0.004755565048, 0.004748141095, 0.004740227096, 0.004731817958, ...
        0.004722909228, 0.004713497184, 0.004703578904, 0.004693152332, ...
        0.004682216325, 0.004670770698, 0.004658816245, 0.004646354755, ...
        0.004633389009, 0.004619922770, 0.004605960752, 0.004591508584, ...
        0.004576572754, 0.004561160551, 0.004545279991, 0.004528939733, ...
        0.004512148992, 0.004494917441, 0.004477255113, 0.004459172292, ...
        0.004440679412, 0.004421786948, 0.004402505313, 0.004382844754, ...
        0.004362815256, 0.004342426454, 0.004321687545, 0.004300607214, ...
        0.004279193573, 0.004257454099, 0.004235395601, 0.004213024178, ...
        0.004190345212, 0.004167363355, 0.004144082539, 0.004120505997, ...
        0.004096636296, 0.004072475378, 0.004048024621, 0.004023284898, ...
        0.003998256657, 0.003972940002, 0.003947334780, 0.003921440679, ...
        0.003895257323, 0.003868784372, 0.003842021626, 0.003814969120, ...
        0.003787627222, 0.003759996728, 0.003732078943, 0.003703875769, ...
        0.003675389765, 0.003646624221, 0.003617583196, 0.003588271569, ...
        0.003558695057, 0.003528860235, 0.003498774536, 0.003468446237, ...
        0.003437884440, 0.003407099027, 0.003376100618, 0.003344900504, ...
        0.003313510581, 0.003281943260, 0.003250211388, 0.003218328139, ...
        0.003186306922, 0.003154161264, 0.003121904707, 0.003089550693, ...
        0.003057112453, 0.003024602902, 0.002992034530, 0.002959419304, ...
        0.002926768579, 0.002894093007, 0.002861402467, 0.002828706000, ...
        0.002796011756, 0.002763326951, 0.002730657846, 0.002698009728, ...
        0.002665386912, 0.002632792754, 0.002600229679, 0.002567699223, ...
        0.002535202083, 0.002502738185, 0.002470306756, 0.002437906415, ...
        0.002405535266, 0.002373190996, 0.002340870987, 0.002308572425, ...
        0.002276292416, 0.002244028101, 0.002211776768, 0.002179535967, ...
        0.002147303617, 0.002115078108, 0.002082858396, 0.002050644089, ...
        0.002018435527, 0.001986233843, 0.001954041023, 0.001921859941, ...
        0.001889694393, 0.001857549108, 0.001825429752, 0.001793342912, ...
        0.001761296070, 0.001729297566, 0.001697356541, 0.001665482879, ...
        0.001633687126, 0.001601980408, 0.001570374340, 0.001538880923, ...
        0.001507512436, 0.001476281334, 0.001445200126, 0.001414281270, ...
        0.001383537057, 0.001352979500, 0.001322620233, 0.001292470407, ...
        0.001262540599, 0.001232840729, 0.001203379980, 0.001174166740, ...
        0.001145208547, 0.001116512045, 0.001088082964, 0.001059926097, ...
        0.001032045304, 0.001004443523, 0.000977122793, 0.000950084291, ...
        0.000923328383, 0.000896854682, 0.000870662116, 0.000844749009, ...
        0.000819113166, 0.000793751965, 0.000768662450, 0.000743841437, ...
        0.000719285612, 0.000694991628, 0.000670956213, 0.000647176257, ...
        0.000623648909, 0.000600371663, 0.000577342430, 0.000554559610, ...
        0.000532022152, 0.000509729602, 0.000487682136, 0.000465880589, ...
        0.000444326464, 0.000423021935, 0.000401969827, 0.000381173598, ...
        0.000360637298, 0.000340365516, 0.000320363327, 0.000300636218, ...
        0.000281190009, 0.000262030769, 0.000243164726, 0.000224598164, ...
        0.000206337332, 0.000188388336, 0.000170757041, 0.000153448970, ...
        0.000136469211, 0.000119822322, 0.000103512249, 0.000087542246, ...
        0.000071914809, 0.000056631617, 0.000041693485, 0.000027100327, ...
        0.000012851138, -0.000001056021, -0.000014624018, -0.000027856630, ...
        -0.000040758515, -0.000053335171, -0.000065592878, -0.000077538630, ...
        -0.000089180064, -0.000100525363, -0.000111583170, -0.000122362479, ...
        -0.000132872532, -0.000143122702, -0.000153122383, -0.000162880871, ...
        -0.000172407255, -0.000181710300, -0.000190798343, -0.000199679191, ...
        -0.000208360028, -0.000216847326, -0.000225146777, -0.000233263221, ...
        -0.000241200601, -0.000248961920, -0.000256549216, -0.000263963554, ...
        -0.000271205024, -0.000278272762, -0.000285164980, -0.000291879010, ...
        -0.000298411360, -0.000304757787, -0.000310913376, -0.000316872631, ...
        -0.000322629576, -0.000328177859, -0.000333510867, -0.000338621843, ...
        -0.000343504007, -0.000348150670, -0.000352555363, -0.000356711946, ...
        -0.000360614724, -0.000364258554, -0.000367638943, -0.000370752142, ...
        -0.000373595219, -0.000376166139, -0.000378463812, -0.000380488141, ...
        -0.000382240056, -0.000383721524, -0.000384935557, -0.000385886201, ...
        -0.000386578507, -0.000387018495, -0.000387213101, -0.000387170118, ...
        -0.000386898114, -0.000386406356, -0.000385704709, -0.000384803542, ...
        -0.000383713619, -0.000382445990, -0.000381011877, -0.000379422563, ...
        -0.000377689276, -0.000375823080, -0.000373834773, -0.000371734777, ...
        -0.000369533056, -0.000367239018, -0.000364861452, -0.000362408450, ...
        -0.000359887364, -0.000357304758, -0.000354666377, -0.000351977139, ...
        -0.000349241120, -0.000346461575, -0.000343640950, -0.000340780925, ...
        -0.000337882450, -0.000334945811, -0.000331970689, -0.000328956235, ...
        -0.000325901155, -0.000322803792, -0.000319662220, -0.000316474338, ...
        -0.000313237965, -0.000309950933, -0.000306611182, -0.000303216849, ...
        -0.000299766348, -0.000296258456, -0.000292692375, -0.000289067796, ...
        -0.000285384955, -0.000281644665, -0.000277848354, -0.000273998076, ...
        -0.000270096521, -0.000266147004, -0.000262153450, -0.000258120361, ...
        -0.000254052771, -0.000249956200, -0.000245836584, -0.000241700207, ...
        -0.000237553622, -0.000233403563, -0.000229256857, -0.000225120328, ...
        -0.000221000704, -0.000216904517, -0.000212838013, -0.000208807055, ...
        -0.000204817042, -0.000200872823, -0.000196978626, -0.000193137994, ...
        -0.000189353730, -0.000185627852, -0.000181961566, -0.000178355242, ...
        -0.000174808415, -0.000171319789, -0.000167887260, -0.000164507953, ...
        -0.000161178271, -0.000157893954, -0.000154650156, -0.000151441530, ...
        -0.000148262318, -0.000145106460, -0.000141967702, -0.000138839714, ...
        -0.000135716209, -0.000132591065, -0.000129458453, -0.000126312954, ...
        -0.000123149681, -0.000119964392, -0.000116753598, -0.000113514662, ...
        -0.000110245887, -0.000106946596, -0.000103617193, -0.000100259219, ...
        -0.000096875386, -0.000093469597, -0.000090046951, -0.000086613733, ...
        -0.000083177384, -0.000079746457, -0.000076330559, -0.000072940270, ...
        -0.000069587057, -0.000066283168, -0.000063041519, -0.000059875563, ...
        -0.000056799157, -0.000053826416, -0.000050971567, -0.000048248791, ...
        -0.000045672075, -0.000043255049, -0.000041010841, -0.000038951924, ...
        -0.000037089977, -0.000035435750, -0.000033998941, -0.000032788081, ...
        -0.000031810440, -0.000031071936, -0.000030577072, -0.000030328878];
tic;
total_input_f_pcm=aec('input_data/rec_mic_0_0_short.pcm', 'input_data/rec_mic_1_0_short.pcm', 'input_data/rec_spk_l_0_short.pcm', filter_coeff);
toc;
disp(['����ʱ��: ',num2str(toc)]);
tmp1 = total_input_f_pcm(:,1:129);
tmp2 = total_input_f_pcm(:,130:258);
[out1, out2] = compose(tmp1, tmp2, filter_coeff);
x = size(out1, 1);
y = size(out1, 2);
out1 = reshape(out1.', [x*y, 1]);
out2 = reshape(out2.', [x*y, 1]);

file_id=fopen('out_aec_0_0.pcm','wb');
fwrite(file_id, out1,'int16');
fclose(file_id);

file_id=fopen('out_aec_1_0.pcm','wb');
fwrite(file_id, out2,'int16');
fclose(file_id);

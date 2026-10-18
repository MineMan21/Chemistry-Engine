// Periodic table data for all 118 elements.
// Columns: Z | symbol | name | atomic mass (u) | Pauling electronegativity (blank = unknown,
// estimated below) | covalent radius (pm) | allowed valences (first = normal, higher ones need a
// more electronegative partner) | CPK/Jmol colour | melting point (K) | boiling point (K)
(function (CE) {
  'use strict';

  const RAW = `
1|H|Hydrogen|1.008|2.20|31|1|FFFFFF|14|20
2|He|Helium|4.0026||28|0|D9FFFF||4
3|Li|Lithium|6.94|0.98|128|1|CC80FF|454|1603
4|Be|Beryllium|9.0122|1.57|96|2|C2FF00|1560|2742
5|B|Boron|10.81|2.04|84|3|FFB5B5|2349|4200
6|C|Carbon|12.011|2.55|76|4|909090|3915|3915
7|N|Nitrogen|14.007|3.04|71|3|3050F8|63|77
8|O|Oxygen|15.999|3.44|66|2|FF0D0D|54|90
9|F|Fluorine|18.998|3.98|57|1|90E050|54|85
10|Ne|Neon|20.180||58|0|B3E3F5|25|27
11|Na|Sodium|22.990|0.93|166|1|AB5CF2|371|1156
12|Mg|Magnesium|24.305|1.31|141|2|8AFF00|923|1363
13|Al|Aluminium|26.982|1.61|121|3|BFA6A6|934|2792
14|Si|Silicon|28.085|1.90|111|4|F0C8A0|1687|3538
15|P|Phosphorus|30.974|2.19|107|3/5|FF8000|317|554
16|S|Sulfur|32.06|2.58|105|2/4/6|FFFF30|388|718
17|Cl|Chlorine|35.45|3.16|102|1/3/5/7|1FF01F|172|239
18|Ar|Argon|39.948||106|0|80D1E3|84|87
19|K|Potassium|39.098|0.82|203|1|8F40D4|337|1032
20|Ca|Calcium|40.078|1.00|176|2|3DFF00|1115|1757
21|Sc|Scandium|44.956|1.36|170|3|E6E6E6|1814|3109
22|Ti|Titanium|47.867|1.54|160|2/3/4|BFC2C7|1941|3560
23|V|Vanadium|50.942|1.63|153|3/4/5|A6A6AB|2183|3680
24|Cr|Chromium|51.996|1.66|139|3/6|8A99C7|2180|2944
25|Mn|Manganese|54.938|1.55|139|2/4/7|9C7AC7|1519|2334
26|Fe|Iron|55.845|1.83|132|2/3|E06633|1811|3134
27|Co|Cobalt|58.933|1.88|126|2/3|F090A0|1768|3200
28|Ni|Nickel|58.693|1.91|124|2|50D050|1728|3186
29|Cu|Copper|63.546|1.90|132|1/2|C88033|1358|2835
30|Zn|Zinc|65.38|1.65|122|2|7D80B0|693|1180
31|Ga|Gallium|69.723|1.81|122|3|C28F8F|303|2673
32|Ge|Germanium|72.630|2.01|120|4|668F8F|1211|3106
33|As|Arsenic|74.922|2.18|119|3/5|BD80E3|887|887
34|Se|Selenium|78.971|2.55|120|2/4/6|FFA100|494|958
35|Br|Bromine|79.904|2.96|120|1/3/5|A62929|266|332
36|Kr|Krypton|83.798|3.00|116|0/2|5CB8D1|116|120
37|Rb|Rubidium|85.468|0.82|220|1|702EB0|312|961
38|Sr|Strontium|87.62|0.95|195|2|00FF00|1050|1655
39|Y|Yttrium|88.906|1.22|190|3|94FFFF|1799|3609
40|Zr|Zirconium|91.224|1.33|175|4|94E0E0|2128|4682
41|Nb|Niobium|92.906|1.6|164|3/5|73C2C9|2750|5017
42|Mo|Molybdenum|95.95|2.16|154|4/6|54B5B5|2896|4912
43|Tc|Technetium|98|1.9|147|4/7|3B9E9E|2430|4538
44|Ru|Ruthenium|101.07|2.2|146|3/4|248F8F|2607|4423
45|Rh|Rhodium|102.91|2.28|142|3|0A7D8C|2237|3968
46|Pd|Palladium|106.42|2.20|139|2/4|006985|1828|3236
47|Ag|Silver|107.87|1.93|145|1|C0C0C0|1235|2435
48|Cd|Cadmium|112.41|1.69|144|2|FFD98F|594|1040
49|In|Indium|114.82|1.78|142|3|A67573|430|2345
50|Sn|Tin|118.71|1.96|139|2/4|668080|505|2875
51|Sb|Antimony|121.76|2.05|139|3/5|9E63B5|904|1860
52|Te|Tellurium|127.60|2.1|138|2/4/6|D47A00|723|1261
53|I|Iodine|126.90|2.66|139|1/3/5/7|940094|387|457
54|Xe|Xenon|131.29|2.6|140|0/2/4/6|429EB0|161|165
55|Cs|Caesium|132.91|0.79|244|1|57178F|302|944
56|Ba|Barium|137.33|0.89|215|2|00C900|1000|2170
57|La|Lanthanum|138.91|1.10|207|3|70D4FF|1193|3737
58|Ce|Cerium|140.12|1.12|204|3/4|FFFFC7|1068|3716
59|Pr|Praseodymium|140.91|1.13|203|3|D9FFC7|1208|3793
60|Nd|Neodymium|144.24|1.14|201|3|C7FFC7|1297|3347
61|Pm|Promethium|145|1.13|199|3|A3FFC7|1315|3273
62|Sm|Samarium|150.36|1.17|198|2/3|8FFFC7|1345|2067
63|Eu|Europium|151.96|1.2|198|2/3|61FFC7|1099|1802
64|Gd|Gadolinium|157.25|1.20|196|3|45FFC7|1585|3546
65|Tb|Terbium|158.93|1.2|194|3|30FFC7|1629|3503
66|Dy|Dysprosium|162.50|1.22|192|3|1FFFC7|1680|2840
67|Ho|Holmium|164.93|1.23|192|3|00FF9C|1734|2993
68|Er|Erbium|167.26|1.24|189|3|00E675|1802|3141
69|Tm|Thulium|168.93|1.25|190|3|00D452|1818|2223
70|Yb|Ytterbium|173.05|1.1|187|2/3|00BF38|1097|1469
71|Lu|Lutetium|174.97|1.27|187|3|00AB24|1925|3675
72|Hf|Hafnium|178.49|1.3|175|4|4DC2FF|2506|4876
73|Ta|Tantalum|180.95|1.5|170|5|4DA6FF|3290|5731
74|W|Tungsten|183.84|2.36|162|4/6|2194D6|3695|6203
75|Re|Rhenium|186.21|1.9|151|4/7|267DAB|3459|5869
76|Os|Osmium|190.23|2.2|144|4/8|266696|3306|5285
77|Ir|Iridium|192.22|2.20|141|3/4|175487|2719|4701
78|Pt|Platinum|195.08|2.28|136|2/4|D0D0E0|2041|4098
79|Au|Gold|196.97|2.54|136|1/3|FFD123|1337|3129
80|Hg|Mercury|200.59|2.00|132|1/2|B8B8D0|234|630
81|Tl|Thallium|204.38|1.62|145|1/3|A6544D|577|1746
82|Pb|Lead|207.2|1.87|146|2/4|575961|601|2022
83|Bi|Bismuth|208.98|2.02|148|3|9E4FB5|545|1837
84|Po|Polonium|209|2.0|140|2/4|AB5C00|527|1235
85|At|Astatine|210|2.2|150|1|754F45|575|610
86|Rn|Radon|222|2.2|150|0/2|428296|202|212
87|Fr|Francium|223|0.79|260|1|420066|300|950
88|Ra|Radium|226|0.9|221|2|007D00|973|2010
89|Ac|Actinium|227|1.1|215|3|70ABFA|1323|3471
90|Th|Thorium|232.04|1.3|206|4|00BAFF|2023|5061
91|Pa|Protactinium|231.04|1.5|200|5|00A1FF|1841|4300
92|U|Uranium|238.03|1.38|196|4/6|008FFF|1405|4404
93|Np|Neptunium|237|1.36|190|5|0080FF|917|4273
94|Pu|Plutonium|244|1.28|187|4|006BFF|913|3501
95|Am|Americium|243|1.13|180|3|545CF2|1449|2880
96|Cm|Curium|247|1.28|169|3|785CE3|1613|3383
97|Bk|Berkelium|247|1.3|168|3|8A4FE3|1259|2900
98|Cf|Californium|251|1.3|168|3|A136D4|1173|1743
99|Es|Einsteinium|252|1.3|165|3|B31FD4|1133|1269
100|Fm|Fermium|257|1.3|167|3|B31FBA|1800|
101|Md|Mendelevium|258|1.3|173|3|B30DA6|1100|
102|No|Nobelium|259|1.3|176|2|BD0D87|1100|
103|Lr|Lawrencium|266|1.3|161|3|C70066|1900|
104|Rf|Rutherfordium|267||157|4|CC0059|2400|5800
105|Db|Dubnium|268||149|5|D1004F||
106|Sg|Seaborgium|269||143|6|D90045||
107|Bh|Bohrium|270||141|7|E00038||
108|Hs|Hassium|277||134|8|E6002E||
109|Mt|Meitnerium|278||129|3|EB0026||
110|Ds|Darmstadtium|281||128|2|ED1A3A||
111|Rg|Roentgenium|282||121|1|EE3354||
112|Cn|Copernicium|285||122|2|F04D6B|283|340
113|Nh|Nihonium|286||136|1|F26682|700|1430
114|Fl|Flerovium|289||143|2|F28099|284|
115|Mc|Moscovium|290||162|1/3|F599B0|670|1400
116|Lv|Livermorium|293||175|2|F5B3C7|709|1085
117|Ts|Tennessine|294||165|1|F7CCDB|723|883
118|Og|Oganesson|294||157|0|F7E0EB|325|450
`;

  const NOBLE = new Set([2, 10, 18, 36, 54, 86, 118]);
  const ALKALI = new Set([3, 11, 19, 37, 55, 87]);
  const ALKALINE = new Set([4, 12, 20, 38, 56, 88]);
  const HALOGEN = new Set([9, 17, 35, 53, 85, 117]);
  const NONMETAL = new Set([1, 6, 7, 8, 15, 16, 34]);
  const METALLOID = new Set([5, 14, 32, 33, 51, 52]);
  const POST_TRANSITION = new Set([13, 31, 49, 50, 81, 82, 83, 84, 113, 114, 115, 116]);

  function periodOf(z) {
    if (z <= 2) return 1;
    if (z <= 10) return 2;
    if (z <= 18) return 3;
    if (z <= 36) return 4;
    if (z <= 54) return 5;
    if (z <= 86) return 6;
    return 7;
  }
  const PERIOD_START = [0, 1, 3, 11, 19, 37, 55, 87];

  // Returns { group, fIndex } where fIndex is set for lanthanides/actinides.
  function groupOf(z) {
    const p = periodOf(z);
    if (p === 1) return { group: z === 1 ? 1 : 18, fIndex: -1 };
    const o = z - PERIOD_START[p] + 1;
    if (p <= 3) return { group: o <= 2 ? o : o + 10, fIndex: -1 };
    if (p <= 5) return { group: o, fIndex: -1 };
    if (o <= 2) return { group: o, fIndex: -1 };
    if (o <= 17) return { group: 3, fIndex: o - 3 };
    return { group: o - 14, fIndex: -1 };
  }

  function categoryOf(z, fIndex) {
    if (NOBLE.has(z)) return 'noble';
    if (ALKALI.has(z)) return 'alkali';
    if (ALKALINE.has(z)) return 'alkaline';
    if (HALOGEN.has(z)) return 'halogen';
    if (NONMETAL.has(z)) return 'nonmetal';
    if (METALLOID.has(z)) return 'metalloid';
    if (POST_TRANSITION.has(z)) return 'post-transition';
    if (fIndex >= 0) return z < 89 ? 'lanthanide' : 'actinide';
    return 'transition';
  }

  const CATEGORY_LABEL = {
    noble: 'Noble gas', alkali: 'Alkali metal', alkaline: 'Alkaline earth metal',
    halogen: 'Halogen', nonmetal: 'Reactive nonmetal', metalloid: 'Metalloid',
    'post-transition': 'Post-transition metal', transition: 'Transition metal',
    lanthanide: 'Lanthanide', actinide: 'Actinide'
  };

  // Lennard-Jones-like well depth (K) for non-metal neutral atoms (roughly UFF-derived).
  const VDW_EPS = {
    1: 22, 2: 4.5, 5: 90, 6: 53, 7: 36, 8: 31, 9: 26, 10: 9.3, 14: 200, 15: 153, 16: 138,
    17: 120, 18: 41, 32: 190, 33: 155, 34: 180, 35: 220, 36: 49, 51: 225, 52: 230,
    53: 300, 54: 74, 85: 330, 86: 55, 118: 67
  };

  const list = [null];
  const bySymbol = {};
  RAW.trim().split('\n').forEach(function (line) {
    const c = line.split('|');
    const z = +c[0];
    const gp = groupOf(z);
    const category = categoryOf(z, gp.fIndex);
    const valences = c[6].split('/').map(Number);
    const el = {
      z: z,
      symbol: c[1],
      name: c[2],
      mass: +c[3],
      en: c[4] === '' ? null : +c[4],
      enEstimated: false,
      covRadius: +c[5],
      valences: valences,
      maxValence: Math.max.apply(null, valences),
      normalValence: valences[0],
      color: '#' + c[7],
      mp: c[8] === '' ? null : +c[8],
      bp: c[9] === '' ? null : +c[9],
      period: periodOf(z),
      group: gp.group,
      fIndex: gp.fIndex,
      category: category,
      categoryLabel: CATEGORY_LABEL[category],
      isMetal: ['alkali', 'alkaline', 'transition', 'post-transition', 'lanthanide', 'actinide'].indexOf(category) >= 0,
      isNoble: category === 'noble',
      radioactive: z === 43 || z === 61 || z >= 84
    };
    // Valence electrons (main group only; used for lone pairs / VSEPR).
    if (el.group <= 2) el.valenceElectrons = el.group;
    else if (el.group >= 13) el.valenceElectrons = el.group - 10;
    else el.valenceElectrons = null;
    if (z === 2) el.valenceElectrons = 2;
    list.push(el);
    bySymbol[el.symbol] = el;
  });

  // Unknown electronegativities (super-heavy elements): estimate from the element one period above.
  for (let z = 1; z <= 118; z++) {
    const el = list[z];
    if (el.en !== null) continue;
    if (el.isNoble && z < 36) { el.en = 0; continue; }
    const above = list[z - 32];
    el.en = above && above.en ? above.en : 1.5;
    el.enEstimated = true;
  }

  // Simulation radius (in world units, 1 unit = 10 pm): a van-der-Waals-like contact radius.
  for (let z = 1; z <= 118; z++) {
    const el = list[z];
    el.radius = Math.min(26, (el.covRadius * 0.75 + 80) / 10);
    el.vdwEps = VDW_EPS[z] || (el.isMetal ? 60 : 120);
  }

  // Elements that naturally exist as diatomic molecules.
  const DIATOMIC = new Set([1, 7, 8, 9, 17, 35, 53, 85, 117]); // (astatine and tennessine as the heavier halogens)

  CE.elements = list;
  CE.bySymbol = bySymbol;
  CE.DIATOMIC = DIATOMIC;
  CE.CATEGORY_LABEL = CATEGORY_LABEL;
})(typeof window !== 'undefined' ? (window.CE = window.CE || {}) : (globalThis.CE = globalThis.CE || {}));

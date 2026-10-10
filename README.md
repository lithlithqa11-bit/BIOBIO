# 🧬 Bio-Impact Analyzer
> **منصة النمذجة الجزيئية المتقدمة ومحاكاة الارتباط الدوائي (Molecular Modeling & Mutational Docking)**

---

## 📌 نظرة عامة (Overview)
منصة **Bio-Impact Analyzer** هي أداة حاسوبية وبيولوجية متكاملة مصممة للصيادلة، أطباء الأورام، والباحثين في علم الأحياء الجزيئي وتصميم الأدوية، تتيح:
1. **التحليل الهيكلي ومقارنة الطفرات**: مقارنة البنية البلورية للبروتين السليم (Wild-Type) مع البروتين المصاب (Mutant)، مع حساب التغير في مساحة السطح المعرضة للمذيب (SASA)، محاذاة التسلسل (Sequence Alignment)، وتصنيف الأثر الكيميائي والفيزيائي للطفرة مع عرض ثلاثي أبعاد تفاعلي بواسطة `py3Dmol`.
2. **محاكاة الإرساء الجزيئي المزدوج (Matched WT-vs-Variant Docking)**: إرساء الجزيئات الدوائية على الجيوب المستهدفة ومقارنة السليم بالطافر مع التحقق الصارم من تطابق التسلسل ($\ge 90\%$) ومحاذاة البقايا المتجانسة وحساب $\Delta \text{Score}$.
3. **الإرساء المتقاطع الاستكشاري (Exploratory Cross-Protein Docking)**: إرساء ومقارنة أي بروتينين مستقلين بشكل متوازٍ ومنفصل دون فرض تجانس جيني أو حساب $\Delta \text{Score}$ مضلل.
4. **الموثوقية والتحقق الإحصائي (Scientific Rigor & Reproducibility)**: دعم البذور العشوائية المتعددة (Multi-seed Docking)، تجميع الوضعيات في عناقيد إحصائية (Pose Clustering)، والتحقق من إعادة الإرساء (Redocking Validation) وتدقيق هوية الربائط عبر قاموس RCSB.

---

## 🏗️ هيكلية المشروع (Project Architecture)

```
BIOBIO/
│
├── app.py                       # نقطة الانطلاق الرئيسية للتطبيق (Streamlit)
│
├── analysis/                    # حزمة التحليل الهيكلي والمقارنة الجزيئية
│   ├── __init__.py              # تصدير دوال الحزمة
│   ├── constants.py             # ثوابت الأحماض الأمينية والخواص الفيزيوكيميائية
│   ├── fetch.py                 # جلب هياكل PDB (Deposited & Assembly) وقاعدة الطفرات
│   ├── structure.py             # قراءة البنية، فحص السلاسل، خوارزميات KDTree، وحسابات SASA
│   ├── alignment.py             # خوارزميات المحاذاة (Biopython) وتحليل أثر الطفرات
│   ├── visualization.py         # العرض ثلاثي الأبعاد (py3Dmol) ورسم Plotly البياني
│   └── ui.py                    # واجهة المستخدم الخاصة بالتحليل الهيكلي وعزل الجلسات
│
├── docking/                     # حزمة الإرساء الجزيئي والتحليل الدوائي
│   ├── __init__.py
│   ├── ui.py                    # واجهة المستخدم للإرساء (Matched & Cross-docking)
│   ├── engine.py                # إدارة تشغيل محرك AutoDock Vina
│   ├── binding_site.py          # حساب صندوق الجيب (Grid Box)
│   ├── receptor_prep.py         # تهيئة وتجهيز البروتين المستقبل وتصحيح أعمدة العناصر
│   ├── ligand_prep.py           # تهيئة الليجاند وتحويل SMILES إلى PDBQT
│   ├── openbabel_io.py          # معالجة آمنة لملفات وتنسيقات OpenBabel
│   ├── pair_validation.py       # التحقق البيولوجي من تطابق الأزواج ومحاذاة بقايا الجيب
│   ├── comparison.py            # المقارنة المزدوجة المتطابقة والإرساء المتقاطع
│   ├── interactions.py          # رصد الروابط الهيدروجينية ومسافات التلامس
│   ├── presets.py               # مكتبة الحالات السريرية البشرية
│   ├── pubchem.py               # جلب الأدوية والمعلومات الجزيئية من PubChem
│   ├── rcsb.py                  # تنزيل ملفات الـ PDB وجلب صيغ الربائط من قاموس RCSB
│   ├── validation.py            # فحص RMSD وتجميع العناقيد والتحقق من هوية الربائط
│   ├── storage.py               # إدارة مسارات التشغيل والمخرجات وتدقيق الـ Manifests
│   └── models.py                # نماذج وهياكل البيانات (Dataclasses)
│
├── data/                        # البيانات والمراجع المحلية
│   └── mutation.json            # قاعدة بيانات ربط البروتينات السليمة بالمصابة تلقائياً
│
├── tools/                       # الأدوات التنفيذية
│   └── vina/
│       ├── vina                 # محرك AutoDock Vina (Linux)
│       └── vina.exe             # محرك AutoDock Vina (Windows)
│
├── tests/                       # حزمة الاختبارات الآلية والتحقق العلمي
│   ├── __init__.py
│   └── test_repairs.py          # اختبارات الوحدة الشاملة (51 اختباراً)
│
├── docking_runs/                # مجلد حفظ مخرجات عمليات الإرساء
│   └── .gitkeep
│
├── requirements.txt             # اعتمادات ومكتبات بايثون
└── .gitignore                   # استثناءات Git
```

---

## 🚀 طريقة التشغيل (How to Run)

### 1. تثبيت الحزم والمتطلبات:
```bash
pip install -r requirements.txt
```

### 2. تشغيل التطبيق:
```bash
streamlit run app.py
```

---

## 🧪 التحقق والاختبار (Verification)

### 1. اختبار استيراد الوحدات الأساسية:
```bash
python -c "import analysis; import docking; print('All modules loaded successfully!')"
```

### 2. تشغيل حزمة الاختبارات الآلية (51 اختباراً علمياً وهيكلياً):
```bash
python -m unittest discover tests
```

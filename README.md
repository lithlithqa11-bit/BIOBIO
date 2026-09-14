# 🧬 Bio-Impact Analyzer
> **منصة النمذجة الجزيئية المتقدمة ومحاكاة الارتباط الدوائي (Molecular Modeling & Mutational Docking)**

---

## 📌 نظرة عامة (Overview)
منصة **Bio-Impact Analyzer** هي أداة حاسوبية وبيولوجية متكاملة مصممة للصيادلة، أطباء الأورام، والباحثين في علم الأحياء الجزيئي وتصميم الأدوية، تتيح:
1. **التحليل الهيكلي ومقارنة الطفرات**: مقارنة البنية البلورية للبروتين السليم (Wild-Type) مع البروتين المصاب (Mutant)، مع حساب التغير في مساحة السطح المعرضة للمذيب (SASA)، محاذاة التسلسل (Sequence Alignment)، وتصنيف الأثر الكيميائي والفيزيائي للطفرة مع عرض ثلاثي أبعاد تفاعلي بواسطة `py3Dmol`.
2. **محاكاة الإرساء الجزيئي (Molecular Docking)**: إرساء الجزيئات الدوائية على الجيوب المستهدفة باستخدام محرك `AutoDock Vina`، مع دعم البذور العشوائية المتعددة (Multi-seed Docking)، تحليل الروابط الهيدروجينية والتلامس، وحساب تقارير الفعالية السريرية ومقاومة الطفرات.

---

## 🏗️ هيكلية المشروع (Project Architecture)

```
Bio-main/
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
│   └── ui.py                    # واجهة المستخدم الخاصة بالتحليل الهيكلي
│
├── docking/                     # حزمة الإرساء الجزيئي والتحليل الدوائي
│   ├── __init__.py
│   ├── ui.py                    # واجهة المستخدم للإرساء
│   ├── engine.py                # إدارة تشغيل محرك AutoDock Vina
│   ├── binding_site.py          # حساب صندوق الجيب (Grid Box)
│   ├── receptor_prep.py         # تهيئة وتجهيز البروتين المستقبل
│   ├── ligand_prep.py           # تهيئة الليجاند وتحويل SMILES إلى PDBQT
│   ├── comparison.py            # المقارنة المزدوجة بين السليم والمصاب
│   ├── interactions.py          # رصد الروابط الهيدروجينية ومسافات التلامس
│   ├── presets.py               # مكتبة الحالات السريرية البشرية
│   ├── pubchem.py               # جلب الأدوية والمعلومات الجزيئية من PubChem
│   ├── rcsb.py                  # تنزيل ملفات الـ PDB
│   ├── validation.py            # التحقق الإحصائي وفحص RMSD
│   ├── storage.py               # إدارة مسارات التشغيل والمخرجات
│   └── models.py                # نماذج وهياكل البيانات (Dataclasses)
│
├── data/                        # البيانات والمراجع المحلية
│   └── mutation.json            # قاعدة بيانات ربط البروتينات السليمة بالمصابة تلقائياً
│
├── tools/                       # الأدوات التنفيذية
│   └── vina/
│       └── vina.exe             # محرك AutoDock Vina
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
لاختبار سلامة كافة موديولات وحزم المشروع:
```bash
python -c "import analysis; import docking; print('All modules loaded successfully!')"
```

# 3D molecular visualization using py3Dmol and SASA comparative Plotly charts
import py3Dmol
import plotly.graph_objects as go
import pandas as pd

def render_protein_3d(pdb_string, bg_color='#111', style_type='cartoon',
                      show_surface=True, surface_opacity=0.3, surface_type='MS',
                      mutations=None, mut_color='red', zoom_to_mutations=False,
                      focus_mut=None, show_ligands=True,
                      ligand_style='أعواد ملونة مع كرات (Sticks & Spheres)'):
    """توليد كود HTML لعرض بنية البروتين ثلاثية الأبعاد باستخدام مكتبة py3Dmol."""
    view = py3Dmol.view(width="100%", height=450)
    view.addModel(pdb_string, 'pdb')
    view.setBackgroundColor(bg_color)

    # تحديد نمط العرض (كرتوني، عصي، كرات)
    style_dict = {style_type: {'color': 'spectrum'}}
    view.setStyle({'model': -1}, style_dict)

    # إظهار جزيئات الليجاند والربائط (HETATM)
    if show_ligands and ligand_style != "إخفاء":
        if "أعواد فقط" in ligand_style or ligand_style == "sticks":
            het_style = {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.35}}
        elif "كرات" in ligand_style and "أعواد" not in ligand_style:
            het_style = {'sphere': {'colorscheme': 'greenCarbon', 'radius': 1.0}}
        else:
            het_style = {
                'stick': {'colorscheme': 'greenCarbon', 'radius': 0.35},
                'sphere': {'colorscheme': 'greenCarbon', 'radius': 0.6}
            }
        view.addStyle({'hetflag': True}, het_style)
    else:
        view.setStyle({'hetflag': True}, {})

    # إظهار السطح الخارجي للبروتين
    if show_surface:
        surf_kind = py3Dmol.MS if surface_type == 'MS' else (py3Dmol.VDW if surface_type == 'VDW' else py3Dmol.SAS)
        view.addSurface(surf_kind, {'opacity': surface_opacity, 'color': '#FFC107'}, {'model': -1})

    # تلوين وتمييز أماكن الطفرات (دمج الأنماط في قاموس واحد لمنع الكتابة الفوقية)
    if mutations:
        for mut in mutations:
            view.addStyle(mut, {
                style_type: {'color': mut_color},
                'stick': {'colorscheme': 'yellowCarbon', 'radius': 0.3},
                'sphere': {'color': mut_color, 'radius': 1.2}
            })

    # التحكم في تقريب الكاميرا (Zoom)
    if focus_mut:
        view.zoomTo(focus_mut)
    elif zoom_to_mutations and mutations:
        view.zoomTo({'or': mutations})
    else:
        view.zoomTo()

    return view._make_html()

def build_sasa_figure(comparison_df: pd.DataFrame) -> go.Figure:
    """بناء رسم بياني تفاعلي باستخدام Plotly لمقارنة منحنيات مساحة SASA وتمييز الطفرات."""
    fig = go.Figure()

    # رسم منحنى البروتين السليم
    fig.add_trace(go.Scatter(
        x=comparison_df['plot_index'],
        y=comparison_df['SASA_H'],
        name='البروتين السليم (WT)',
        line=dict(color='#00ff88', width=2),
        fill=None,
        customdata=comparison_df[['res_num_h', 'السليم']].values,
        hovertemplate="🟢 السليم: %{customdata[1]} (رقم %{customdata[0]}) — SASA: %{y:.2f} Å²<extra></extra>"
    ))

    # رسم منحنى البروتين المصاب مع تظليل الفرق
    fig.add_trace(go.Scatter(
        x=comparison_df['plot_index'],
        y=comparison_df['SASA_M'],
        name='البروتين المصاب (MT)',
        line=dict(color='#ff3333', width=2),
        fill='tonexty', 
        fillcolor='rgba(255, 51, 51, 0.1)',
        customdata=comparison_df[['res_num_m', 'المصاب']].values,
        hovertemplate="🔴 المصاب: %{customdata[1]} (رقم %{customdata[0]}) — SASA: %{y:.2f} Å²<extra></extra>"
    ))

    # تمييز مواقع الطفرات بنجوم صفراء على الرسم
    mutations_df = comparison_df[comparison_df['السليم'] != comparison_df['المصاب']]
    if not mutations_df.empty:
        fig.add_trace(go.Scatter(
            x=mutations_df['plot_index'],
            y=mutations_df['SASA_M'],
            mode='markers',
            name='مواقع الطفرات',
            marker=dict(color='yellow', size=10, symbol='star', line=dict(color='black', width=1)),
            hovertemplate="⭐ طفرة: %{customdata[2]} (%{customdata[0]}) ← %{customdata[3]} (%{customdata[1]})<br>التأثير: %{customdata[4]}<extra></extra>",
            customdata=mutations_df[['res_num_h', 'res_num_m', 'السليم', 'المصاب', 'Impact']].values
        ))

    # وضع أرقام الأحماض الحقيقية كـ Labels أسفل الرسم
    tick_step = max(1, len(comparison_df) // 20)
    fig.update_layout(
        template="plotly_dark",
        height=450,
        hovermode="x unified",
        xaxis=dict(
            title="موقع الحمض الأميني (Residue Position)",
            rangeslider=dict(visible=True),
            tickmode='array',
            tickvals=comparison_df['plot_index'][::tick_step].tolist(),
            ticktext=comparison_df['res_num_label'][::tick_step].astype(str).tolist()
        ),
        yaxis=dict(title="SASA (Å²)"),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        margin=dict(l=0, r=0, t=30, b=0)
    )
    return fig

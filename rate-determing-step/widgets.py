""" Programmes pour créer des widgets

Une class Widget (et FloatSlider, IntSlider) qui permet
d'enregistrer des paramètres avant de créer le widget. La 
syntaxe est la même que les ipywidgets. Ceci permet entre
autre de mettre de séparer le fond (les parmètres) de la forme
(l'axe), ce que matplotlib ne permet pas.

Des fonctions : 
    make_param_widgets, 
    make_choose_plot, 
    make_reset_button, 
    make_log_button,
    make_start_stop_animation

"""

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, CheckButtons

class Widget(object):
    value = None
    description = ""
    def __init__(self, **kwd):
        for key, val in kwd.items():
            if not hasattr(type(self), key):
                raise Exception('Attribut "{}" non valide'.format(key))
            setattr(self, key, val)

class FloatSlider(Widget):
    min = 0
    max = 1

    def make_mpl_widget(self, ax, update):
        w = Slider(ax, self.description, self.min, self.max, valinit=self.value)
        w.on_changed(update)
        return w

class IntSlider(Widget):
    min = 1
    max = 10

    def make_mpl_widget(self, ax, update):
        w = Slider(ax, self.description, self.min, self.max, valinit=self.value, valstep=1)
        w.on_changed(update)
        return w


#class Checkbox(Widget):
#    def make_mpl_widget(self, ax, update):
#        pos = list(ax.get_position().bounds)
#        pos[2] = 0.1
#        ax.set_position(pos)
#        w = CheckButtons(ax, [self.description], actives=[self.value])
#        w.on_clicked(update)
#        return w


slider_color = 'lightgoldenrodyellow'


def make_param_widgets(parameters, plot_data, slider_box):
    """ Crée automatiquement les widget matplotlib

    parameters : dictionnaire contenant les parameters
    plot_data : callback fonction
    slider_box : boite dans laquelle mettre les sliders
    """
    f = plt.gcf()
    n = len(parameters)
    x0, y0, W, H = slider_box
    height = H/n

    def update(val=None):
        values = {}
        for key, w in mpl_widgets.items():
            values[key] = w.val
        plot_data(**values)
        plt.draw()


#    default = {key:val.value for key, val in parameters.items()}
    mpl_widgets = {}
    for i, (key, elm) in enumerate(parameters.items()):
        ax = f.add_axes([x0, y0+height*i, W, height], facecolor=slider_color)
        mpl_widgets[key] = elm.make_mpl_widget(ax, update)        
    update()

    return mpl_widgets
    

def make_choose_plot(lines, box, which=None):
    """Create check button for the lines

    lines : all the lines
    box : the box where the check button will be placed
    which : by default (None) : all the plots
        otherwise : key of the lines 
            or tuple of keys (one button for many lines)
    """
    f = plt.gcf()
    ax = f.add_axes(box, facecolor=slider_color)

    if which is None:
        which = lines.keys()
    labels = []
    is_active = []
    for key in which:
        if isinstance(key, tuple):
            key = key[0]
        elm = lines[key]
        labels.append(key)
        is_active.append(elm.get_visible())
    choose_widget = CheckButtons(ax, labels, is_active)

    def chooseplot(label):
        states = choose_widget.get_status() 
        for i, key in enumerate(which):
            state = states[i]
            if isinstance(key, tuple):
                for k in key:
                    lines[k].set_visible(state)
            else:
                lines[key].set_visible(state)            
        f.canvas.draw_idle()

    choose_widget.on_clicked(chooseplot)

    return choose_widget

def make_reset_button(mpl_widgets, box=[0.8, 0.005, 0.1, 0.04]):
    f = plt.gcf()
    ax = f.add_axes(box, facecolor=slider_color)

    button = Button(ax, 'Reset', color=slider_color, hovercolor='0.975')

    def reset(event):
        for widget in mpl_widgets.values():
            widget.reset() 

    button.on_clicked(reset) # Lorsqu'on clique sur "reset", on applique la fonction reset definie au dessus

    return button

def make_log_button(ax, box=[0.015, 0.05, 0.12, 0.15], ylims=None):
    """ Make a log button

    ax : the axis
    ylims : None or dictionary with the ylim for 'linear' and 'log' scales
    """
    f = plt.gcf()
    ax_btn = f.add_axes(box, facecolor=slider_color)

    labels = ['log x', 'log y']
    widget = CheckButtons(ax_btn, labels, [ax.get_xscale()=='log', ax.get_yscale()=='log'])

    def set_log(label):
        if label =='log x':
            method = 'set_xscale'
            index = 0
        if label =='log y':
            method = 'set_yscale'
            index = 1
        state = 'log' if widget.get_status()[index] else 'linear'
        getattr(ax, method)(state)
        if ylims is not None:
            if label=='log y':
                ax.set_ylim(ylims[state])

        f.canvas.draw_idle()

    widget.on_clicked(set_log)

    return widget

def make_start_stop_animation(anim, box=[0.015, 0.05, 0.12, 0.1], start_animation=True):
    f = plt.gcf()
    ax_btn = f.add_axes(box, facecolor=slider_color)

    labels = ['anim']
    widget = CheckButtons(ax_btn, labels, [start_animation])

    def set_anim(label):
        if widget.get_status()[0]:
            print('on')
            anim.event_source.start()
        else:
            print('off')
            anim.event_source.stop()

    widget.on_clicked(set_anim)

    return widget
def justify_paragraph(string, width=40):
    out = ['']
    in_equation = False
    for word in string.split(' '):
        if word.startswith('$'):
            in_equation = True
        if not in_equation:
            if len(out[-1]) + len(word) + 1> width:
                out.append('')
        out[-1] += ' ' + word
        if word.endswith('$'):
            in_equation = False
    out = [elm for elm in out if elm.strip()] # Remove blanck lines
    return '\n'.join(out)

def justify(string, width=40):
    paragraphs = ['']
    for lines in string.split('\n'):
        if lines.strip()=='':
            paragraphs.append('')
        else:
            paragraphs[-1] += lines.strip() + ' '

    out = [justify_paragraph(paragraph, width=width) for paragraph in paragraphs]
    return '\n\n'.join(out)
    


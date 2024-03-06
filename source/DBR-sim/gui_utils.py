"""Helper functions and -classes for propeller GUI."""

from PySide6 import QtCore, QtWidgets, QtGui


KEYCODES = {
    "Enter": 16777220,
    "Escape": 16777216,
    "UpArrow": 16777235,
    "DownArrow": 16777237,
}


class Container(QtWidgets.QWidget):
    """Container with layout to which given child widgets are added."""

    _layout = QtWidgets.QLayout

    def __init__(self, parent=None, children=[]):
        super(Container, self).__init__(parent=parent)

        self.layout = self._layout(self)
        for child in children:
            if isinstance(child, QtWidgets.QWidget):
                self.layout.addWidget(child)
            elif isinstance(child, QtWidgets.QSpacerItem):
                self.layout.addItem(child)
            elif child == "stretch":
                self.layout.addStretch()


class HLayout(QtWidgets.QHBoxLayout):
    """HBoxLayout without contents margins."""

    def __init__(self, parent=None):
        super(HLayout, self).__init__(parent)
        self.setContentsMargins(0, 0, 0, 0)


class VLayout(QtWidgets.QVBoxLayout):
    """VBoxLayout without contents margins."""

    def __init__(self, parent=None):
        super(VLayout, self).__init__(parent)
        self.setContentsMargins(0, 0, 0, 0)


class HContainer(Container):
    """Container with HLayout to which given child widgets are added."""

    _layout = HLayout

    def __init__(self, parent=None, children=[]):
        super(HContainer, self).__init__(parent=parent, children=children)


class VContainer(Container):
    """Container with VLayout to which given child widgets are added."""

    _layout = VLayout

    def __init__(self, parent=None, children=[]):
        super(VContainer, self).__init__(parent=parent, children=children)


class ListSettingWidget(QtWidgets.QWidget):
    """Listwidget connected to 'Add' and 'Remove' buttons and an input field.

    This class can be used as a base to inherit from and add customized functionality.
    Overwrite the class variables to change:
    - New item field widget
    - New item field placeholder text
    - Text in the 'Add' button
    - Text in the 'Remove' button
    """

    _NEW_ITEM_FIELD = QtWidgets.QLineEdit
    _PLACEHOLDER_TEXT = ""
    _ADD_BTN_TEXT = "Add"
    _REMOVE_BTN_TEXT = "Remove"
    def _IS_ITEM_NAME_VALID(*args, name=None): True

    def __init__(self, parent=None):
        super(ListSettingWidget, self).__init__(parent=parent)

        layout = VLayout(self)

        self.new_item_field = self._NEW_ITEM_FIELD()
        self.new_item_field.setPlaceholderText(self._PLACEHOLDER_TEXT)

        self.add_btn = QtWidgets.QPushButton(self._ADD_BTN_TEXT)
        self.add_btn.setFixedWidth(60)
        header = HContainer(children=[self.new_item_field, self.add_btn])

        self.listwidget = QtWidgets.QListWidget()
        self.remove_btn = QtWidgets.QPushButton(self._REMOVE_BTN_TEXT)
        self.remove_btn.setFixedWidth(60)
        removal_container = VContainer(children=[self.remove_btn, "stretch"])
        list_container = HContainer(children=[self.listwidget, removal_container])

        layout.addWidget(header)
        layout.addWidget(list_container)

        self.remove_btn.clicked.connect(self._on_remove_clicked)
        self.add_btn.clicked.connect(self._on_add_clicked)
        self.new_item_field.textChanged.connect(self._on_text_changed)
        self.listwidget.itemClicked.connect(self._on_selection_changed)

        self.refresh()

    def _on_selection_changed(self, *args):
        self.refresh()

    def _on_add_clicked(self, *args):
        self.add_item()
        self.new_item_field.setText("")
        self.new_item_field.setPlaceholderText(self._PLACEHOLDER_TEXT)
        self.refresh()

    def _on_text_changed(self, *args):
        self.refresh()

    def _on_remove_clicked(self, *args):
        self.remove_selected_items()
        self.refresh()

    def add_item(self):
        item_name = self.new_item_field.text()
        validation_result = self._IS_ITEM_NAME_VALID(name=item_name)
        if validation_result is True:
            self.listwidget.insertItem(0, item_name)
        elif isinstance(validation_result, str):
            print(validation_result)
            tooltip = QtWidgets.QToolTip(validation_result)
            position = self.position(self.new_item_field)
            tooltip.showText(position, validation_result)

    def position(self, widget):
        return widget.mapToGlobal(QtCore.QPoint())

    def remove_selected_items(self):
        for item in self.listwidget.selectedItems():
            row = self.listwidget.row(item)
            self.listwidget.takeItem(row)

    def refresh(self):
        if self.new_item_field.text():
            self.add_btn.setEnabled(True)
        else:
            self.add_btn.setEnabled(False)

        selected_items = self.listwidget.selectedItems()
        if selected_items:
            self.remove_btn.setEnabled(True)
        else:
            self.remove_btn.setEnabled(False)

        self.new_item_field.setFocus()

    def keyPressEvent(self, event):
        enter_pressed = event.key() == KEYCODES.get("Enter")
        if enter_pressed:
            self._on_add_clicked()

        modifiers = event.modifiers()
        shift_pressed = QtCore.Qt.ShiftModifier & modifiers    # TODO: Use these events 
        control_pressed = QtCore.Qt.ShiftModifier & modifiers  # to support selection of
                                                               # multiple items


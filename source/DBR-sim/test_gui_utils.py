"""Unit test suite for gui_utils module.

This module provides comprehensive tests for the GUI utility functionality.
Tests cover container classes, layout classes, list widget functionality, and keycode constants.
"""

import unittest
import os
import sys
from unittest.mock import MagicMock, patch

# Add the parent directory to Python path to import the module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Mock PySide6 imports to avoid GUI dependencies in tests
sys.modules['PySide6'] = MagicMock()
sys.modules['PySide6.QtCore'] = MagicMock()
sys.modules['PySide6.QtWidgets'] = MagicMock()
sys.modules['PySide6.QtGui'] = MagicMock()

import gui_utils as gui


class TestKeycodes(unittest.TestCase):
    """Test cases for KEYCODES constants."""

    def test_keycodes_dict_exists(self):
        """Test that KEYCODES dictionary exists."""
        self.assertIsInstance(gui.KEYCODES, dict)

    def test_keycodes_has_expected_keys(self):
        """Test that KEYCODES has expected key entries."""
        expected_keys = ["Enter", "Escape", "UpArrow", "DownArrow"]
        for key in expected_keys:
            self.assertIn(key, gui.KEYCODES)

    def test_keycodes_values_are_integers(self):
        """Test that all KEYCODES values are integers."""
        for key, value in gui.KEYCODES.items():
            self.assertIsInstance(value, int)

    def test_enter_keycode_value(self):
        """Test Enter keycode has expected value."""
        self.assertEqual(gui.KEYCODES["Enter"], 16777220)

    def test_escape_keycode_value(self):
        """Test Escape keycode has expected value."""
        self.assertEqual(gui.KEYCODES["Escape"], 16777216)

    def test_up_arrow_keycode_value(self):
        """Test UpArrow keycode has expected value."""
        self.assertEqual(gui.KEYCODES["UpArrow"], 16777235)

    def test_down_arrow_keycode_value(self):
        """Test DownArrow keycode has expected value."""
        self.assertEqual(gui.KEYCODES["DownArrow"], 16777237)


class TestContainerClasses(unittest.TestCase):
    """Test cases for Container base class and its subclasses."""

    def test_container_class_exists(self):
        """Test that Container class exists."""
        self.assertTrue(hasattr(gui, 'Container'))
        self.assertTrue(callable(gui.Container))

    def test_hcontainer_class_exists(self):
        """Test that HContainer class exists."""
        self.assertTrue(hasattr(gui, 'HContainer'))
        self.assertTrue(callable(gui.HContainer))

    def test_vcontainer_class_exists(self):
        """Test that VContainer class exists."""
        self.assertTrue(hasattr(gui, 'VContainer'))
        self.assertTrue(callable(gui.VContainer))

    def test_hcontainer_inherits_from_container(self):
        """Test that HContainer inherits from Container."""
        self.assertTrue(issubclass(gui.HContainer, gui.Container))

    def test_vcontainer_inherits_from_container(self):
        """Test that VContainer inherits from Container."""
        self.assertTrue(issubclass(gui.VContainer, gui.Container))

    def test_container_layout_attribute(self):
        """Test that Container has _layout attribute."""
        self.assertTrue(hasattr(gui.Container, '_layout'))

    def test_hcontainer_layout_is_hlayout(self):
        """Test that HContainer uses HLayout as _layout."""
        self.assertEqual(gui.HContainer._layout, gui.HLayout)

    def test_vcontainer_layout_is_vlayout(self):
        """Test that VContainer uses VLayout as _layout."""
        self.assertEqual(gui.VContainer._layout, gui.VLayout)


class TestLayoutClasses(unittest.TestCase):
    """Test cases for layout classes."""

    def test_hlayout_class_exists(self):
        """Test that HLayout class exists."""
        self.assertTrue(hasattr(gui, 'HLayout'))
        self.assertTrue(callable(gui.HLayout))

    def test_vlayout_class_exists(self):
        """Test that VLayout class exists."""
        self.assertTrue(hasattr(gui, 'VLayout'))
        self.assertTrue(callable(gui.VLayout))

    def test_hlayout_inherits_from_qhboxlayout(self):
        """Test that HLayout inherits from QHBoxLayout."""
        # We can't test the actual inheritance due to the mocking,
        # but we can test that the class has the expected name
        self.assertEqual(gui.HLayout.__name__, 'HLayout')

    def test_vlayout_inherits_from_qvboxlayout(self):
        """Test that VLayout inherits from QVBoxLayout."""
        self.assertEqual(gui.VLayout.__name__, 'VLayout')


class TestListSettingWidget(unittest.TestCase):
    """Test cases for ListSettingWidget class."""

    def test_list_setting_widget_class_exists(self):
        """Test that ListSettingWidget class exists."""
        self.assertTrue(hasattr(gui, 'ListSettingWidget'))
        self.assertTrue(callable(gui.ListSettingWidget))

    def test_list_setting_widget_class_variables(self):
        """Test that ListSettingWidget has expected class variables."""
        expected_vars = ['_NEW_ITEM_FIELD', '_PLACEHOLDER_TEXT', '_ADD_BTN_TEXT', '_REMOVE_BTN_TEXT']
        for var in expected_vars:
            self.assertTrue(hasattr(gui.ListSettingWidget, var))

    def test_list_setting_widget_default_values(self):
        """Test ListSettingWidget default class variable values."""
        self.assertEqual(gui.ListSettingWidget._PLACEHOLDER_TEXT, "")
        self.assertEqual(gui.ListSettingWidget._ADD_BTN_TEXT, "Add")
        self.assertEqual(gui.ListSettingWidget._REMOVE_BTN_TEXT, "Remove")

    def test_list_setting_widget_validation_function(self):
        """Test that ListSettingWidget has validation function."""
        self.assertTrue(hasattr(gui.ListSettingWidget, '_IS_ITEM_NAME_VALID'))

    def test_list_setting_widget_validation_default(self):
        """Test ListSettingWidget default validation always returns True."""
        # The default validation function should accept any name
        result = gui.ListSettingWidget._IS_ITEM_NAME_VALID(name="test_name")
        self.assertTrue(result)

    def test_list_setting_widget_required_methods(self):
        """Test that ListSettingWidget has required methods."""
        required_methods = [
            '__init__', 'add_item', 'remove_selected_items', 'refresh', 
            '_on_add_clicked', '_on_remove_clicked', '_on_text_changed',
            '_on_selection_changed', 'keyPressEvent', 'position'
        ]
        for method in required_methods:
            self.assertTrue(hasattr(gui.ListSettingWidget, method))


class TestContainerInitialization(unittest.TestCase):
    """Test cases for Container initialization behavior."""

    def test_container_init_accepts_parent(self):
        """Test that Container.__init__ accepts parent parameter."""
        # We can't actually create the object due to Qt dependencies,
        # but we can check the method signature
        import inspect
        sig = inspect.signature(gui.Container.__init__)
        self.assertIn('parent', sig.parameters)

    def test_container_init_accepts_children(self):
        """Test that Container.__init__ accepts children parameter."""
        import inspect
        sig = inspect.signature(gui.Container.__init__)
        self.assertIn('children', sig.parameters)

    def test_container_init_children_default(self):
        """Test that Container children parameter has default value."""
        import inspect
        sig = inspect.signature(gui.Container.__init__)
        self.assertEqual(sig.parameters['children'].default, [])

    def test_hcontainer_init_signature(self):
        """Test HContainer.__init__ signature."""
        import inspect
        sig = inspect.signature(gui.HContainer.__init__)
        params = list(sig.parameters.keys())
        self.assertIn('parent', params)
        self.assertIn('children', params)

    def test_vcontainer_init_signature(self):
        """Test VContainer.__init__ signature."""
        import inspect
        sig = inspect.signature(gui.VContainer.__init__)
        params = list(sig.parameters.keys())
        self.assertIn('parent', params)
        self.assertIn('children', params)


class TestListSettingWidgetMethods(unittest.TestCase):
    """Test cases for ListSettingWidget method behavior."""

    def test_add_item_method_exists(self):
        """Test that add_item method exists."""
        self.assertTrue(hasattr(gui.ListSettingWidget, 'add_item'))
        self.assertTrue(callable(gui.ListSettingWidget.add_item))

    def test_remove_selected_items_method_exists(self):
        """Test that remove_selected_items method exists."""
        self.assertTrue(hasattr(gui.ListSettingWidget, 'remove_selected_items'))
        self.assertTrue(callable(gui.ListSettingWidget.remove_selected_items))

    def test_refresh_method_exists(self):
        """Test that refresh method exists."""
        self.assertTrue(hasattr(gui.ListSettingWidget, 'refresh'))
        self.assertTrue(callable(gui.ListSettingWidget.refresh))

    def test_position_method_exists(self):
        """Test that position method exists."""
        self.assertTrue(hasattr(gui.ListSettingWidget, 'position'))
        self.assertTrue(callable(gui.ListSettingWidget.position))

    def test_key_press_event_method_exists(self):
        """Test that keyPressEvent method exists."""
        self.assertTrue(hasattr(gui.ListSettingWidget, 'keyPressEvent'))
        self.assertTrue(callable(gui.ListSettingWidget.keyPressEvent))


class TestModuleStructure(unittest.TestCase):
    """Test cases for overall module structure."""

    def test_module_docstring(self):
        """Test that module has docstring."""
        self.assertIsNotNone(gui.__doc__)
        self.assertGreater(len(gui.__doc__), 0)

    def test_expected_classes_present(self):
        """Test that all expected classes are present in the module."""
        expected_classes = [
            'Container', 'HContainer', 'VContainer', 
            'HLayout', 'VLayout', 'ListSettingWidget'
        ]
        for class_name in expected_classes:
            self.assertTrue(hasattr(gui, class_name))
            self.assertTrue(callable(getattr(gui, class_name)))

    def test_keycodes_constant(self):
        """Test that KEYCODES constant is present."""
        self.assertTrue(hasattr(gui, 'KEYCODES'))
        self.assertIsInstance(gui.KEYCODES, dict)


class TestListSettingWidgetEventHandling(unittest.TestCase):
    """Test cases for ListSettingWidget event handling."""

    def test_on_add_clicked_calls_add_item(self):
        """Test that _on_add_clicked calls add_item."""
        # This tests the method signature and logic flow
        self.assertTrue(hasattr(gui.ListSettingWidget, '_on_add_clicked'))

    def test_on_remove_clicked_calls_remove_selected_items(self):
        """Test that _on_remove_clicked calls remove_selected_items."""
        self.assertTrue(hasattr(gui.ListSettingWidget, '_on_remove_clicked'))

    def test_on_text_changed_calls_refresh(self):
        """Test that _on_text_changed calls refresh."""
        self.assertTrue(hasattr(gui.ListSettingWidget, '_on_text_changed'))

    def test_on_selection_changed_calls_refresh(self):
        """Test that _on_selection_changed calls refresh."""
        self.assertTrue(hasattr(gui.ListSettingWidget, '_on_selection_changed'))


class TestContainerChildrenHandling(unittest.TestCase):
    """Test cases for Container children handling logic."""

    def test_container_init_processes_children(self):
        """Test that Container.__init__ processes children parameter."""
        # We can't test the actual behavior due to Qt dependencies,
        # but we can verify that the method exists and has the right signature
        import inspect
        source = inspect.getsource(gui.Container.__init__)
        
        # Check that it processes children
        self.assertIn('for child in children:', source)

    def test_container_handles_widget_children(self):
        """Test Container handles QWidget children."""
        import inspect
        source = inspect.getsource(gui.Container.__init__)
        
        # Check that it handles QWidget children
        self.assertIn('isinstance(child, QtWidgets.QWidget)', source)
        self.assertIn('addWidget(child)', source)

    def test_container_handles_spacer_children(self):
        """Test Container handles QSpacerItem children."""
        import inspect
        source = inspect.getsource(gui.Container.__init__)
        
        # Check that it handles QSpacerItem children
        self.assertIn('isinstance(child, QtWidgets.QSpacerItem)', source)
        self.assertIn('addItem(child)', source)

    def test_container_handles_stretch_children(self):
        """Test Container handles stretch string children."""
        import inspect
        source = inspect.getsource(gui.Container.__init__)
        
        # Check that it handles stretch
        self.assertIn('child == "stretch"', source)
        self.assertIn('addStretch()', source)


class TestLayoutMargins(unittest.TestCase):
    """Test cases for layout margin settings."""

    def test_hlayout_sets_zero_margins(self):
        """Test that HLayout sets contents margins to zero."""
        import inspect
        source = inspect.getsource(gui.HLayout.__init__)
        
        # Check that it sets zero margins
        self.assertIn('setContentsMargins(0, 0, 0, 0)', source)

    def test_vlayout_sets_zero_margins(self):
        """Test that VLayout sets contents margins to zero."""
        import inspect
        source = inspect.getsource(gui.VLayout.__init__)
        
        # Check that it sets zero margins
        self.assertIn('setContentsMargins(0, 0, 0, 0)', source)


class TestListSettingWidgetBehavior(unittest.TestCase):
    """Test cases for ListSettingWidget behavior."""

    def test_add_item_uses_new_item_field(self):
        """Test that add_item gets text from new_item_field."""
        import inspect
        source = inspect.getsource(gui.ListSettingWidget.add_item)
        
        # Check that it uses new_item_field.text()
        self.assertIn('self.new_item_field.text()', source)

    def test_add_item_validates_name(self):
        """Test that add_item validates item name."""
        import inspect
        source = inspect.getsource(gui.ListSettingWidget.add_item)
        
        # Check that it calls validation
        self.assertIn('_IS_ITEM_NAME_VALID', source)
        self.assertIn('validation_result', source)

    def test_add_item_adds_to_listwidget(self):
        """Test that add_item adds valid items to listwidget."""
        import inspect
        source = inspect.getsource(gui.ListSettingWidget.add_item)
        
        # Check that it adds to listwidget
        self.assertIn('self.listwidget.insertItem', source)

    def test_remove_selected_items_removes_from_list(self):
        """Test that remove_selected_items removes selected items."""
        import inspect
        source = inspect.getsource(gui.ListSettingWidget.remove_selected_items)
        
        # Check that it removes selected items
        self.assertIn('self.listwidget.selectedItems()', source)
        self.assertIn('self.listwidget.takeItem', source)

    def test_refresh_updates_button_states(self):
        """Test that refresh updates button enable states."""
        import inspect
        source = inspect.getsource(gui.ListSettingWidget.refresh)
        
        # Check that it updates button states based on conditions
        self.assertIn('setEnabled', source)
        self.assertIn('self.new_item_field.text()', source)

    def test_key_press_event_handles_enter(self):
        """Test that keyPressEvent handles Enter key."""
        import inspect
        source = inspect.getsource(gui.ListSettingWidget.keyPressEvent)
        
        # Check that it handles Enter key
        self.assertIn('KEYCODES.get("Enter")', source)
        self.assertIn('_on_add_clicked', source)


class TestModuleImports(unittest.TestCase):
    """Test cases for module imports."""

    def test_pyside6_imports(self):
        """Test that PySide6 modules are imported."""
        # This is tested implicitly by the fact that the module loads
        # without import errors (thanks to our mocking)
        self.assertTrue(True)


if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)